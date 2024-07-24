/*
        The MIT License (MIT)

        Copyright (c) 2021 José Antonio Fernández Fernández

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.
*/
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace tmd
{
/* ==========================================  DECLARATIONS  ========================================== */
// Small vector 3D class class
template <typename FLOAT>
class Vec3r
{
  public:
    std::array<FLOAT, 3> v;

    Vec3r(){};
    template <typename FLOAT_I>
    Vec3r(const FLOAT_I &x, const FLOAT_I &y, const FLOAT_I &z)
    {
        v[0] = static_cast<FLOAT>(x);
        v[1] = static_cast<FLOAT>(y);
        v[2] = static_cast<FLOAT>(z);
    }
    template <typename SIZE_T>
    const FLOAT &operator[](const SIZE_T &i) const { return v[i]; }
    template <typename SIZE_T>
    FLOAT &operator[](const SIZE_T &i) { return v[i]; }
    FLOAT dot(const Vec3r &u) const { return v[0] * u[0] + v[1] * u[1] + v[2] * u[2]; }
    Vec3r<FLOAT> cross(const Vec3r &u) const { return Vec3r(v[1] * u[2] - v[2] * u[1], -v[0] * u[2] + v[2] * u[0], v[0] * u[1] - v[1] * u[0]); }
    Vec3r<FLOAT> operator+(const Vec3r &u) const { return Vec3r(v[0] + u[0], v[1] + u[1], v[2] + u[2]); }
    Vec3r<FLOAT> operator-(const Vec3r &u) const { return Vec3r(v[0] - u[0], v[1] - u[1], v[2] - u[2]); }
    void operator+=(const Vec3r &u)
    {
        v[0] += u[0];
        v[1] += u[1];
        v[2] += u[2];
    }
    template <typename FLOAT_I>
    Vec3r<FLOAT> operator*(const FLOAT_I &a) const { return Vec3r(static_cast<FLOAT>(a) * v[0], static_cast<FLOAT>(a) * v[1], static_cast<FLOAT>(a) * v[2]); }
    template <typename FLOAT_I>
    Vec3r<FLOAT> operator/(const FLOAT_I &a) const { return Vec3r(v[0] / static_cast<FLOAT>(a), v[1] / static_cast<FLOAT>(a), v[2] / static_cast<FLOAT>(a)); }
    template <typename FLOAT_I>
    void operator/=(const FLOAT_I &a)
    {
        v[0] /= static_cast<FLOAT>(a);
        v[1] /= static_cast<FLOAT>(a);
        v[2] /= static_cast<FLOAT>(a);
    }
    FLOAT squaredNorm() const { return this->dot(*this); }
    FLOAT norm() const { return std::sqrt(this->squaredNorm()); }
    Vec3r<FLOAT> normalized() const { return (*this) / this->norm(); }
    void normalize()
    {
        const FLOAT norm = this->norm();
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }
};
template <typename FLOAT, typename FLOAT_I>
static inline Vec3r<FLOAT> operator*(const FLOAT_I &a, const Vec3r<FLOAT> &v) { return Vec3r<FLOAT>(static_cast<FLOAT>(a) * v[0], static_cast<FLOAT>(a) * v[1], static_cast<FLOAT>(a) * v[2]); }
using Vec3d = Vec3r<double>;
// -----------------------------------

// Point-Triangle distance declarations
enum class NearestEntity
{
    V0,
    V1,
    V2,
    E01,
    E12,
    E02,
    F
};

/**
 * Computes the squared distance, the nearest entity (vertex, edge or face) and the nearest point from a point to a triangle.
 */
static void point_triangle_sq_unsigned(double &distance_sq, NearestEntity &nearest_entity, Vec3d &barycentric, Vec3d &nearest_point, const Vec3d &point, const Vec3d &v0, const Vec3d &v1, const Vec3d &v2);
// -----------------------------------

// Struct that contains the result of a distance query
struct Result
{
    double distance = std::numeric_limits<double>::max();
    Vec3d nearest_point;
    tmd::NearestEntity nearest_entity;
    int triangle_id = -1;
    Vec3d barycentric;
};
// -----------------------------------

/**
 * A class to compute signed and unsigned distances to a connected
 * and watertight triangle mesh.
 */
class TriangleMeshDistance
{
  private:
    /* Private declarations */
    struct BoundingSphere
    {
        Vec3d center;
        double radius;
    };

    struct Node
    {
        BoundingSphere bv_left;
        BoundingSphere bv_right;
        int left = -1; // If left == -1, right is the triangle_id
        int right = -1;
    };

    struct Triangle
    {
        std::array<Vec3d, 3> vertices;
        int id = -1;
    };

    /* Private fields */
    std::vector<Vec3d> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<Node> nodes;
    std::vector<Vec3d> pseudonormals_triangles;
    std::vector<std::array<Vec3d, 3>> pseudonormals_edges;
    std::vector<Vec3d> pseudonormals_vertices;
    BoundingSphere root_bv;
    bool is_constructed = false;
    bool is_manifold = false;

    /* Private methods */
    void _construct();
    void _build_tree(const int node_id, BoundingSphere &bounding_sphere, std::vector<Triangle> &triangles, const int begin, const int end);
    void _query(Result &result, const Node &node, const Vec3d &point) const;

  public:
    /* Public methods */
    TriangleMeshDistance() = default;

    /**
     * @brief Constructs a new TriangleMeshDistance object.
     *
     * @param vertices Pointer to the vertices coordinates array in xyzxyz... layout.
     * @param n_vertices Number of vertices.
     * @param triangles Pointer to the conectivity array in ijkijk... layout.
     * @param n_triangles Number of triangles.
     */
    template <typename FLOAT, typename INT, typename SIZE_T>
    TriangleMeshDistance(const FLOAT *vertices, const SIZE_T n_vertices, const INT *triangles, const SIZE_T n_triangles);

    /**
     * @brief Constructs a new TriangleMeshDistance object.
     *
     * @param vertices Vertices of the triangle mesh. Y coordinate of the 3rd vertex should be access by `vertices[2][1]`.
     * @param triangles Triangles of the triangle mesh. Index of the 2nd vertex of the 3rd triangle should be access by `triangles[2][1]`.
     */
    template <typename IndexableVector3double, typename IndexableVector3int>
    TriangleMeshDistance(const std::vector<IndexableVector3double> &vertices, const std::vector<IndexableVector3int> &triangles);

    /**
     * @brief Initializes an existing TriangleMeshDistance object (including empty ones).
     *
     * @param vertices Pointer to the vertices coordinates array in xyzxyz... layout.
     * @param n_vertices Number of vertices.
     * @param triangles Pointer to the conectivity array in ijkijk... layout.
     * @param n_triangles Number of triangles.
     */
    template <typename FLOAT, typename INT, typename SIZE_T>
    void construct(const FLOAT *vertices, const SIZE_T n_vertices, const INT *triangles, const SIZE_T n_triangles);

    /**
     * @brief Initializes an existing TriangleMeshDistance object (including empty ones).
     *
     * @param vertices Vertices of the triangle mesh. Y coordinate of the 3rd vertex should be access by `vertices[2][1]`.
     * @param triangles Triangles of the triangle mesh. Index of the 2nd vertex of the 3rd triangle should be access by `triangles[2][1]`.
     */
    template <typename IndexableVector3double, typename IndexableVector3int>
    void construct(const std::vector<IndexableVector3double> &vertices, const std::vector<IndexableVector3int> &triangles);

    /**
     * @brief Returns the result of the manifold check done at construction time.
     *
     * @return True if the mesh is watertight, with all edges being shared by exactly two triangles, false otherwise.
     */
    bool is_mesh_manifold() const;

    /**
     * @brief Computes the unsigned distance from a point to the triangle mesh. Thread safe.
     *
     * @param point to query from. Typed to `Vec3d` but can be passed as `{x, y, z}`.
     *
     * @return Result containing distance, nearest point on the mesh, nearest entity and the nearest triangle index.
     */
    template <typename IndexableVector3double>
    Result unsigned_distance(const IndexableVector3double &point) const;
    Result unsigned_distance(const std::array<double, 3> &point) const;

    /**
     * @brief Computes the unsigned distance from a point to the triangle mesh. Thread safe.
     *
     * @param point to query from. Typed to `Vec3d` but can be passed as `{x, y, z}`.
     *
     * @return Result containing distance, nearest point on the mesh, nearest entity and the nearest triangle index.
     */
    template <typename IndexableVector3double>
    Result signed_distance(const IndexableVector3double &point) const;
    Result signed_distance(const std::array<double, 3> &point) const;
};
} // namespace tmd

/* ==========================================  DEFINITIONS  ========================================== */
template <typename FLOAT, typename INT, typename SIZE_T>
inline tmd::TriangleMeshDistance::TriangleMeshDistance(const FLOAT *vertices, const SIZE_T n_vertices, const INT *triangles, const SIZE_T n_triangles)
{
    this->construct(vertices, n_vertices, triangles, n_triangles);
}

template <typename IndexableVector3double, typename IndexableVector3int>
inline tmd::TriangleMeshDistance::TriangleMeshDistance(const std::vector<IndexableVector3double> &vertices, const std::vector<IndexableVector3int> &triangles)
{
    this->construct(vertices, triangles);
}

template <typename FLOAT, typename INT, typename SIZE_T>
inline void tmd::TriangleMeshDistance::construct(const FLOAT *vertices, const SIZE_T n_vertices, const INT *triangles, const SIZE_T n_triangles)
{
    this->vertices.resize((size_t)n_vertices);
    for (size_t i = 0; i < (size_t)n_vertices; i++)
    {
        this->vertices[i][0] = (double)vertices[3 * i + 0];
        this->vertices[i][1] = (double)vertices[3 * i + 1];
        this->vertices[i][2] = (double)vertices[3 * i + 2];
    }

    this->triangles.resize((size_t)n_triangles);
    for (size_t i = 0; i < (size_t)n_triangles; i++)
    {
        this->triangles[i][0] = (int)triangles[3 * i + 0];
        this->triangles[i][1] = (int)triangles[3 * i + 1];
        this->triangles[i][2] = (int)triangles[3 * i + 2];
    }
    this->_construct();
}

template <typename IndexableVector3double, typename IndexableVector3int>
inline void tmd::TriangleMeshDistance::construct(const std::vector<IndexableVector3double> &vertices, const std::vector<IndexableVector3int> &triangles)
{
    this->vertices.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++)
    {
        this->vertices[i][0] = (double)vertices[i][0];
        this->vertices[i][1] = (double)vertices[i][1];
        this->vertices[i][2] = (double)vertices[i][2];
    }
    this->triangles.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++)
    {
        this->triangles[i][0] = (int)triangles[i][0];
        this->triangles[i][1] = (int)triangles[i][1];
        this->triangles[i][2] = (int)triangles[i][2];
    }
    this->_construct();
}

inline bool tmd::TriangleMeshDistance::is_mesh_manifold() const
{
    return this->is_manifold;
}

inline tmd::Result tmd::TriangleMeshDistance::signed_distance(const std::array<double, 3> &point) const
{
    const Vec3d p(point[0], point[1], point[2]);
    Result result = this->unsigned_distance(p);

    const std::array<int, 3> &triangle = this->triangles[result.triangle_id];
    Vec3d pseudonormal(0.0, 0.0, 0.0);
    switch (result.nearest_entity)
    {
    case tmd::NearestEntity::V0:
        pseudonormal = this->pseudonormals_vertices[triangle[0]];
        break;
    case tmd::NearestEntity::V1:
        pseudonormal = this->pseudonormals_vertices[triangle[1]];
        break;
    case tmd::NearestEntity::V2:
        pseudonormal = this->pseudonormals_vertices[triangle[2]];
        break;
    case tmd::NearestEntity::E01:
        pseudonormal = this->pseudonormals_edges[result.triangle_id][0];
        break;
    case tmd::NearestEntity::E12:
        pseudonormal = this->pseudonormals_edges[result.triangle_id][1];
        break;
    case tmd::NearestEntity::E02:
        pseudonormal = this->pseudonormals_edges[result.triangle_id][2];
        break;
    case tmd::NearestEntity::F:
        pseudonormal = this->pseudonormals_triangles[result.triangle_id];
        break;

    default:
        break;
    }

    const Vec3d nearest_point(
        result.barycentric[0] * this->vertices[triangle[0]] +
        result.barycentric[1] * this->vertices[triangle[1]] +
        result.barycentric[2] * this->vertices[triangle[2]]);
    const Vec3d u = p - nearest_point;
    result.distance *= (u.dot(pseudonormal) >= 0.0) ? 1.0 : -1.0;

    return result;
}

template <typename IndexableVector3double>
inline tmd::Result tmd::TriangleMeshDistance::signed_distance(const IndexableVector3double &point) const
{
    return this->signed_distance({static_cast<double>(point[0]), static_cast<double>(point[1]), static_cast<double>(point[2])});
}

inline tmd::Result tmd::TriangleMeshDistance::unsigned_distance(const std::array<double, 3> &point) const
{
    if (!this->is_constructed)
    {
        std::cout << "DistanceTriangleMesh error: not constructed." << std::endl;
        exit(-1);
    }

    const Vec3d p(point[0], point[1], point[2]);
    Result result;
    result.distance = std::numeric_limits<double>::max();
    this->_query(result, this->nodes[0], p);
    return result;
}

template <typename IndexableVector3double>
inline tmd::Result tmd::TriangleMeshDistance::unsigned_distance(const IndexableVector3double &point) const
{
    return this->unsigned_distance({static_cast<double>(point[0]), static_cast<double>(point[1]), static_cast<double>(point[2])});
}

inline void tmd::TriangleMeshDistance::_construct()
{
    if (this->triangles.size() == 0)
    {
        std::cout << "DistanceTriangleMesh error: Empty triangle list." << std::endl;
        exit(-1);
    }

    // Build the tree containing the triangles
    std::vector<Triangle> triangles;

    triangles.resize(this->triangles.size());
    for (int i = 0; i < (int)this->triangles.size(); i++)
    {
        triangles[i].id = i;

        const std::array<int, 3> &triangle = this->triangles[i];
        triangles[i].vertices[0] = this->vertices[triangle[0]];
        triangles[i].vertices[1] = this->vertices[triangle[1]];
        triangles[i].vertices[2] = this->vertices[triangle[2]];
    }

    this->nodes.push_back(Node());
    this->_build_tree(0, this->root_bv, triangles, 0, (int)triangles.size());

    // Compute pseudonormals
    //// Edge data structure
    std::unordered_map<uint64_t, Vec3d> edge_normals;
    std::unordered_map<uint64_t, int> edges_count;
    const uint64_t n_vertices = (uint64_t)this->vertices.size();
    auto add_edge_normal = [&](const int i, const int j, const Vec3d &triangle_normal)
    {
        const uint64_t key = std::min(i, j) * n_vertices + std::max(i, j);
        if (edge_normals.find(key) == edge_normals.end())
        {
            edge_normals[key] = triangle_normal;
            edges_count[key] = 1;
        }
        else
        {
            edge_normals[key] += triangle_normal;
            edges_count[key] += 1;
        }
    };
    auto get_edge_normal = [&](const int i, const int j)
    {
        const uint64_t key = std::min(i, j) * n_vertices + std::max(i, j);
        return edge_normals.find(key)->second;
    };

    //// Compute
    this->pseudonormals_triangles.resize(this->triangles.size());
    this->pseudonormals_edges.resize(this->triangles.size());
    this->pseudonormals_vertices.resize(this->vertices.size(), {0, 0, 0});
    for (int i = 0; i < (int)this->triangles.size(); i++)
    {

        // Triangle
        const std::array<int, 3> &triangle = this->triangles[i];
        const Vec3d &a = this->vertices[triangle[0]];
        const Vec3d &b = this->vertices[triangle[1]];
        const Vec3d &c = this->vertices[triangle[2]];

        const Vec3d triangle_normal = (b - a).cross(c - a).normalized();
        this->pseudonormals_triangles[i] = triangle_normal;

        // Vertex
        const double alpha_0 = std::acos((b - a).normalized().dot((c - a).normalized()));
        const double alpha_1 = std::acos((a - b).normalized().dot((c - b).normalized()));
        const double alpha_2 = std::acos((b - c).normalized().dot((a - c).normalized()));
        this->pseudonormals_vertices[triangle[0]] += alpha_0 * triangle_normal;
        this->pseudonormals_vertices[triangle[1]] += alpha_1 * triangle_normal;
        this->pseudonormals_vertices[triangle[2]] += alpha_2 * triangle_normal;

        // Edge
        add_edge_normal(triangle[0], triangle[1], triangle_normal);
        add_edge_normal(triangle[1], triangle[2], triangle_normal);
        add_edge_normal(triangle[0], triangle[2], triangle_normal);
    }

    for (Vec3d &n : this->pseudonormals_vertices)
    {
        n.normalize();
    }

    for (int tri_i = 0; tri_i < (int)this->triangles.size(); tri_i++)
    {
        const std::array<int, 3> &triangle = this->triangles[tri_i];
        this->pseudonormals_edges[tri_i][0] = get_edge_normal(triangle[0], triangle[1]).normalized();
        this->pseudonormals_edges[tri_i][1] = get_edge_normal(triangle[1], triangle[2]).normalized();
        this->pseudonormals_edges[tri_i][2] = get_edge_normal(triangle[0], triangle[2]).normalized();
    }

    // Check that the mesh is watertight: All edges appear exactly twice.
    this->is_manifold = true;
    for (const auto edge_count : edges_count)
    {
        if (edge_count.second != 2)
        {
            this->is_manifold = false;
        }
    }

    this->is_constructed = true;
}

inline void tmd::TriangleMeshDistance::_build_tree(const int node_id, BoundingSphere &bounding_sphere, std::vector<Triangle> &triangles, const int begin, const int end)
{
    const int n_triangles = end - begin;

    if (n_triangles == 0)
    {
        std::cout << "DistanceTriangleMesh::_construct error: Empty leave." << std::endl;
        exit(-1);
    }
    else if (n_triangles == 1)
    {
        // Build node leaf
        this->nodes[node_id].left = -1;
        this->nodes[node_id].right = triangles[begin].id;

        //// Bounding sphere
        const Triangle &tri = triangles[begin];
        const Vec3d center = (tri.vertices[0] + tri.vertices[1] + tri.vertices[2]) / 3.0;
        const double radius = std::max(std::max((tri.vertices[0] - center).norm(), (tri.vertices[1] - center).norm()), (tri.vertices[2] - center).norm());
        bounding_sphere.center = center;
        bounding_sphere.radius = radius;
    }
    else
    {
        // Compute AxisAligned Bounding Box center and largest dimension of all current triangles
        Vec3d top = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
        Vec3d bottom = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        Vec3d center = {0, 0, 0};
        for (int tri_i = begin; tri_i < end; tri_i++)
        {
            for (int vertex_i = 0; vertex_i < 3; vertex_i++)
            {
                const Vec3d &p = triangles[tri_i].vertices[vertex_i];
                center += p;

                for (int coord_i = 0; coord_i < 3; coord_i++)
                {
                    top[coord_i] = std::max(top[coord_i], p[coord_i]);
                    bottom[coord_i] = std::min(bottom[coord_i], p[coord_i]);
                }
            }
        }
        center /= 3 * n_triangles;
        const Vec3d diagonal = top - bottom;
        const int split_dim = (int)(std::max_element(&diagonal[0], &diagonal[0] + 3) - &diagonal[0]);

        // Set node bounding sphere
        double radius_sq = 0.0;
        for (int tri_i = begin; tri_i < end; tri_i++)
        {
            for (int i = 0; i < 3; i++)
            {
                radius_sq = std::max(radius_sq, (center - triangles[tri_i].vertices[i]).squaredNorm());
            }
        }
        bounding_sphere.center = center;
        bounding_sphere.radius = std::sqrt(radius_sq);

        // Sort the triangles according to their center along the split dimension
        std::sort(triangles.begin() + begin, triangles.begin() + end,
                  [split_dim](const Triangle &a, const Triangle &b)
                  {
                      return a.vertices[0][split_dim] < b.vertices[0][split_dim];
                  });

        // Children
        const int mid = (int)(0.5 * (begin + end));

        this->nodes[node_id].left = (int)this->nodes.size();
        this->nodes.push_back(Node());
        this->_build_tree(this->nodes[node_id].left, this->nodes[node_id].bv_left, triangles, begin, mid);

        this->nodes[node_id].right = (int)this->nodes.size();
        this->nodes.push_back(Node());
        this->_build_tree(this->nodes[node_id].right, this->nodes[node_id].bv_right, triangles, mid, end);
    }
}

inline void tmd::TriangleMeshDistance::_query(Result &result, const Node &node, const Vec3d &point) const
{
    // End of recursion
    if (node.left == -1)
    {
        const int triangle_id = node.right;
        const std::array<int, 3> &triangle = this->triangles[node.right]; // If left == -1, right is the triangle_id
        const Vec3d &v0 = this->vertices[triangle[0]];
        const Vec3d &v1 = this->vertices[triangle[1]];
        const Vec3d &v2 = this->vertices[triangle[2]];

        double distance_sq;
        tmd::NearestEntity nearest_entity;
        Vec3d barycentric;
        Vec3d nearest_point;
        tmd::point_triangle_sq_unsigned(distance_sq, nearest_entity, barycentric, nearest_point, point, v0, v1, v2);

        if (distance_sq < result.distance * result.distance)
        {
            result.nearest_entity = nearest_entity;
            result.nearest_point = nearest_point;
            result.barycentric = barycentric;
            result.distance = std::sqrt(distance_sq);
            result.triangle_id = triangle_id;
        }
    }

    // Recursion
    else
    {
        // Find which child bounding volume is closer
        const double d_left = (point - node.bv_left.center).norm() - node.bv_left.radius;
        const double d_right = (point - node.bv_right.center).norm() - node.bv_right.radius;

        if (d_left < d_right)
        {

            // Overlap test
            if (d_left < result.distance)
            {
                this->_query(result, this->nodes[node.left], point);
            }

            if (d_right < result.distance)
            {
                this->_query(result, this->nodes[node.right], point);
            }
        }
        else
        {
            if (d_right < result.distance)
            {
                this->_query(result, this->nodes[node.right], point);
            }
            if (d_left < result.distance)
            {
                this->_query(result, this->nodes[node.left], point);
            }
        }
    }
}

static void tmd::point_triangle_sq_unsigned(double &distance_sq, NearestEntity &nearest_entity, Vec3d &barycentric, Vec3d &nearest_point, const Vec3d &p, const Vec3d &a, const Vec3d &b, const Vec3d &c)
{
    // This function is a modified version of the one found in the Real-Time Collision Detection book by Ericson.
    Vec3d ab = b - a;
    Vec3d ac = c - a;
    Vec3d bc = c - b;

    // Compute parametric position s for projection P’ of P on AB
    double snom = (p - a).dot(ab), sdenom = (p - b).dot(a - b);
    // Compute parametric position t for projection P’ of P on AC
    double tnom = (p - a).dot(ac), tdenom = (p - c).dot(a - c);
    if (snom <= 0.0 && tnom <= 0.0)
    {
        nearest_entity = NearestEntity::V0;
        barycentric = {1.0, 0.0, 0.0};
        nearest_point = a;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }

    // Compute parametric position u for projection P’ of P on BC
    double unom = (p - b).dot(bc), udenom = (p - c).dot(b - c);
    if (sdenom <= 0.0 && unom <= 0.0)
    {
        nearest_entity = NearestEntity::V1;
        barycentric = {0.0, 1.0, 0.0};
        nearest_point = b;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }
    if (tdenom <= 0.0 && udenom <= 0.0)
    {
        nearest_entity = NearestEntity::V2;
        barycentric = {0.0, 0.0, 1.0};
        nearest_point = c;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }

    // Normal for the triangle
    Vec3d n = ab.cross(ac);

    // Check if P is outside AB
    double vc = n.dot((a - p).cross(b - p));
    if (vc <= 0.0 && snom >= 0.0 && sdenom >= 0.0)
    {
        double arc = snom / (snom + sdenom);
        nearest_entity = NearestEntity::E01;
        barycentric = {1.0 - arc, arc, 0.0};
        nearest_point = barycentric[0] * a + barycentric[1] * b;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }

    // Check if P is outside BC
    double va = n.dot((b - p).cross(c - p));
    if (va <= 0.0 && unom >= 0.0 && udenom >= 0.0)
    {
        double arc = unom / (unom + udenom);
        nearest_entity = NearestEntity::E12;
        barycentric = {0.0, 1.0 - arc, arc};
        nearest_point = barycentric[1] * b + barycentric[2] * c;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }

    // Check if P is outside AC
    double vb = n.dot((c - p).cross(a - p));
    if (vb <= 0.0 && tnom >= 0.0 && tdenom >= 0.0)
    {
        double arc = tnom / (tnom + tdenom);
        nearest_entity = NearestEntity::E02;
        barycentric = {1.0 - arc, 0.0, arc};
        nearest_point = barycentric[0] * a + barycentric[2] * c;
        distance_sq = (p - nearest_point).squaredNorm();
        return;
    }

    // P must project inside the triangle; compute using barycentric coordinates
    double u = va / (va + vb + vc);
    double v = vb / (va + vb + vc);
    double w = 1.0 - u - v; // = vc / (va + vb + vc)
    nearest_entity = NearestEntity::F;
    barycentric = {u, v, w};
    nearest_point = u * a + v * b + w * c;
    distance_sq = (p - nearest_point).squaredNorm();
    return;
}
