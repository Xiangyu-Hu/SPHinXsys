# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate example mesh with:
# julia -e 'include("mesh/mesh.jl"); generate_ring_mesh(filename="rings.msh")'

using Gmsh: gmsh
using LinearAlgebra

"""
    generate_ring_mesh(;
        filename::AbstractString,
        wire_width                         = 1.0,
        inner_radius                       = 10.0,
        outer_radius                       = 100.0,
        rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
        rot_θ::Real                        = π / 2,
        verbose::Integer                   = 5,
        gui::Bool                          = false
    )

Generate a mesh for the rings example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - wire_width - width of the rings
  - inner_radius - radius of the inner ring
  - outer_radius - radius of the outer ring
  - rot_center - center of rotation
  - rot_axis - axis of rotation
  - rot_θ - angle of rotation about rot_axis, originating at rot_center
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_ring_mesh(;
    filename::AbstractString,
    wire_width                         = 1.0,
    inner_radius                       = 10.0,
    outer_radius                       = 100.0,
    rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
    rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
    rot_θ::Real                       = π / 6,
    verbose::Integer                   = 5,
    gui::Bool                          = false
)
    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "rings" in gmsh.model.list()
        gmsh.model.setCurrent("rings")
        gmsh.model.remove()
    end
    gmsh.model.add("rings")

    # Geometry parameters (in μm)
    farfield_radius = 10.0 * outer_radius

    # Mesh parameters
    l_ring = 2.0
    l_farfield = 200.0

    # Origin
    p0 = kernel.addPoint(0.0, 0.0, 0.0)

    # Inner ring
    h0 = 0.5 * wire_width
    r1 = inner_radius - h0
    r2 = inner_radius + h0
    x1 = sqrt(r1^2 - h0^2)
    x2 = sqrt(r2^2 - h0^2)

    pi1 = kernel.addPoint(x1, -h0, 0.0, l_ring)
    pi2 = kernel.addPoint(x1, h0, 0.0, l_ring)
    pi3 = kernel.addPoint(x2, h0, 0.0, l_ring)
    pi4 = kernel.addPoint(x2, -h0, 0.0, l_ring)

    li1 = kernel.addLine(pi1, pi2)
    li2 = kernel.addLine(pi2, pi3)
    li3 = kernel.addLine(pi3, pi4)
    li4 = kernel.addLine(pi4, pi1)

    inner_terminal_loop = kernel.addCurveLoop([li1, li2, li3, li4])
    inner_terminal = kernel.addPlaneSurface([inner_terminal_loop])

    pi5 = kernel.addPoint(-r1, 0.0, 0.0, l_ring)
    pi6 = kernel.addPoint(-r2, 0.0, 0.0, l_ring)

    ai1 = kernel.addCircleArc(pi2, p0, pi5)
    ai2 = kernel.addCircleArc(pi5, p0, pi1)
    ai3 = kernel.addCircleArc(pi3, p0, pi6)
    ai4 = kernel.addCircleArc(pi6, p0, pi4)

    inner_ring_loop = kernel.addCurveLoop([ai1, ai2, -li4, -ai4, -ai3, -li2])
    inner_ring = kernel.addPlaneSurface([inner_ring_loop])

    # Outer ring
    r1 = outer_radius - h0
    r2 = outer_radius + h0
    x1 = sqrt(r1^2 - h0^2)
    x2 = sqrt(r2^2 - h0^2)

    po1 = kernel.addPoint(x1, -h0, 0.0, l_ring)
    po2 = kernel.addPoint(x1, h0, 0.0, l_ring)
    po3 = kernel.addPoint(x2, h0, 0.0, l_ring)
    po4 = kernel.addPoint(x2, -h0, 0.0, l_ring)

    lo1 = kernel.addLine(po1, po2)
    lo2 = kernel.addLine(po2, po3)
    lo3 = kernel.addLine(po3, po4)
    lo4 = kernel.addLine(po4, po1)

    outer_terminal_loop = kernel.addCurveLoop([lo1, lo2, lo3, lo4])
    outer_terminal = kernel.addPlaneSurface([outer_terminal_loop])

    po5 = kernel.addPoint(-r1, 0.0, 0.0, l_ring)
    po6 = kernel.addPoint(-r2, 0.0, 0.0, l_ring)

    ao1 = kernel.addCircleArc(po2, p0, po5)
    ao2 = kernel.addCircleArc(po5, p0, po1)
    ao3 = kernel.addCircleArc(po3, p0, po6)
    ao4 = kernel.addCircleArc(po6, p0, po4)

    outer_ring_loop = kernel.addCurveLoop([ao1, ao2, -lo4, -ao4, -ao3, -lo2])
    outer_ring = kernel.addPlaneSurface([outer_ring_loop])

    # Auxiliary surfaces
    inner_gap_loop = kernel.addCurveLoop([ai1, ai2, li1])
    inner_gap = kernel.addPlaneSurface([inner_gap_loop])

    outer_gap_loop_in = kernel.addCurveLoop([ai3, ai4, -li3])
    outer_gap_in = kernel.addPlaneSurface([outer_gap_loop_in])
    outer_gap_loop_out = kernel.addCurveLoop([ao1, ao2, lo1])
    outer_gap_out = kernel.addPlaneSurface([outer_gap_loop_out])
    outer_gap, _ = kernel.cut([(2, outer_gap_out)], [(2, outer_gap_in)])
    @assert length(outer_gap) == 1
    outer_gap = first(outer_gap)[2]

    # Add external box
    domain = kernel.addBox(
        -farfield_radius,
        -farfield_radius,
        -farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius
    )

    # Apply a rotation transformation to all entities in the model
    rot_axis ./= norm(rot_axis)
    kernel.rotate(
        kernel.getEntities(),
        rot_center[1],
        rot_center[2],
        rot_center[3],
        rot_axis[1],
        rot_axis[2],
        rot_axis[3],
        rot_θ
    )

    kernel.synchronize()

    # Add physical groups
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")

    _, farfield_boundaries = gmsh.model.getAdjacencies(3, domain)
    farfield_group = gmsh.model.addPhysicalGroup(2, farfield_boundaries, -1, "farfield")

    rings_group = gmsh.model.addPhysicalGroup(2, [inner_ring, outer_ring], -1, "rings")

    inner_terminal_group =
        gmsh.model.addPhysicalGroup(2, [inner_terminal], -1, "terminal_inner")
    outer_terminal_group =
        gmsh.model.addPhysicalGroup(2, [outer_terminal], -1, "terminal_outer")

    inner_gap_group = gmsh.model.addPhysicalGroup(2, [inner_gap], -1, "hole_inner")
    outer_gap_group = gmsh.model.addPhysicalGroup(2, [outer_gap], -1, "hole_outer")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_ring)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(
        1,
        "SurfacesList",
        [inner_ring, outer_ring, inner_terminal, outer_terminal]
    )
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", 6.0 * outer_radius)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    mesh_curves = last.(
        gmsh.model.getBoundary(
            [(2, x) for x in [inner_ring, outer_ring, inner_terminal, outer_terminal]],
            true,
            false,
            false
        )
    )

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", mesh_curves)
    gmsh.model.mesh.field.setNumber(2, "Sampling", 30)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", l_ring)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "DistMax", 6.0 * outer_radius)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.model.mesh.embed(
        2,
        [inner_terminal, inner_ring, outer_terminal, outer_ring, inner_gap, outer_gap],
        3,
        domain
    )

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundaries: ", farfield_group)
    println("Ring boundaries: ", rings_group)
    println("Inner terminal: ", inner_terminal_group)
    println("Outer terminal: ", outer_terminal_group)
    println("Inner hole: ", inner_gap_group)
    println("Outer hole: ", outer_gap_group)
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
