import bpy

class SPHinXsysGlobalProps(bpy.types.PropertyGroup):
    m_SimDomainLower: bpy.props.FloatVectorProperty(name = 'Simulation Domain Lower', description = 'Simulation Domain Lower', default = (-1, -1, -1))
    m_SimDomainUpper: bpy.props.FloatVectorProperty(name = 'Simulation Domain Upper', description = 'Simulation Domain Upper', default = (1, 1, 1))
    m_ParticleReferenceSpacing: bpy.props.FloatProperty(name = 'Fluid Reference Spacing', description = 'Fluid Reference Spacing', default = 0.05)
    m_Gravity: bpy.props.FloatProperty(name = 'Gravity', description = 'Gravity', default = 1.0)
    m_EndTime: bpy.props.FloatProperty(name = 'End Time', description = 'End Time', default = 20.0)
    m_Rho0f: bpy.props.FloatProperty(name = 'Rho0 Ref', description = 'Rho0 Ref', default = 1.0)
    m_URef: bpy.props.FloatProperty(name = 'Vel Ref', description = 'Vel Ref', default = 20.0)
    m_Cf: bpy.props.FloatProperty(name = 'Sound Ref', description = 'Sound Ref', default = 200.0, )
    m_WaterBlockFilePath: bpy.props.StringProperty(name = 'Water Block File Path', description = 'Water Block File Path', default = "./Mesh/Fluids/")
    m_RigidBlockFilePath: bpy.props.StringProperty(name = 'Rigid Block File Path', description = 'Rigid Block File Path', default = "./Mesh/Rigids/")
    m_ElasticBlockFilePath: bpy.props.StringProperty(name = 'Elastic Block FilePath', description = 'Elastic Block File Path', default = "./Mesh/Elastics/")
    # m_WaterBlockFileName: bpy.props.StringProperty(name = 'Water Block File Name', description = 'Water Block File Name', default = "water.stl")
    # m_RigidBlockFileName: bpy.props.StringProperty(name = 'Rigid Block File Name', description = 'Rigid Block File Name', default = "rigid.stl")
    # m_ElasticBlockFileName: bpy.props.StringProperty(name = 'Elastic Block Name', description = 'Elastic Block File Name', default = "elastic.stl")
    
# items = [('RigidSolid', 'Rigid Solid', 'This object will be defined as a rigid solid.', 0), ('ElasticSolid', 'Elastic Solid', 'This object will be defined as a elastic solid.', 1), ('Fluid', 'Fluid', 'This object will be defined as a fluid.', 2)]
class SPHinXsysProps(bpy.types.PropertyGroup):
    m_TypeofMedia: bpy.props.EnumProperty(name = 'Type of Media', items = 
                                          [('RigidSolid', 'Rigid Solid', 'This object will be defined as a rigid solid.', 0), 
                                           ('ElasticSolid', 'Elastic Solid', 'This object will be defined as a elastic solid.', 1), 
                                           ('Fluid', 'Fluid', 'This object will be defined as a fluid.', 2)], 
                                           description='This property defines the media type of the current mesh.')
    m_Export: bpy.props.BoolProperty(name = 'Export', description = 'Whether the mesh needs to be exported.', default = True)
    # m_FluidBlocksName

    m_FluidDensity: bpy.props.FloatProperty(name = 'Fluid Density', description = 'Fluid Density', default = 1.0)
    m_FluidCharacteristicVelocity: bpy.props.FloatProperty(name = 'Fluid Characteristic Velocity', description = 'Fluid Characteristic Velocity', default = 10.0)
    m_FluidSpeedofSound: bpy.props.FloatProperty(name = 'Fluid Speed of Sound', description = 'Fluid Speed of Sound', default = 10.0)
    m_FluidReynoldsNumber: bpy.props.FloatProperty(name = 'Fluid Reynolds Number', description = 'Fluid Reynolds Number', default = 10.0)
    m_FluidDynamicsViscosity: bpy.props.FloatProperty(name = 'Fluid Dynamics Viscosity', description = 'Fluid Dynamics Viscosity', default = 10.0)

    m_SolidReferenceDensity: bpy.props.FloatProperty(name = 'Solid Reference Density', description = 'Solid Reference Density', default = 10.0)
    m_SolidYoungsModulus: bpy.props.FloatProperty(name = 'Solid Youngs Modulus', description = 'Solid Youngs Modulus', default = 1.4e3)
    m_SolidPoissonRatio: bpy.props.FloatProperty(name = 'Solid Poisson Ratio', description = 'Solid Poisson Ratio', default = 0.4)

    

def registerSPHinXsysProperties():
        bpy.types.Scene.m_SPHinXsysGlobalSettings = bpy.props.PointerProperty(type = SPHinXsysGlobalProps)
        bpy.types.Mesh.m_SPHinXsysSettings = bpy.props.PointerProperty(type = SPHinXsysProps)

m_FluidBlocksName = []
m_RigidBlocksName = []
m_ElasticBlocksName = []

def appendBlocksName(blockType, blockName):
    if(blockType == 'Fluid'):
        if(not(blockName in m_FluidBlocksName)):
            m_FluidBlocksName.append(blockName)
    elif(blockType == 'RigidSolid'):
        if(not(blockName in m_RigidBlocksName)):
            m_RigidBlocksName.append(blockName)
    elif(blockType == 'ElasticSolid'):
        if(not(blockName in m_ElasticBlocksName)):
            m_ElasticBlocksName.append(blockName)

def deleteBlocksName(blockType, blockName):
    if(blockType == 'Fluid'):
        if(blockName in m_FluidBlocksName):
            m_FluidBlocksName.remove(blockName)
    elif(blockType == 'RigidSolid'):
        if(blockName in m_RigidBlocksName):
            m_RigidBlocksName.remove(blockName)
    elif(blockType == 'ElasticSolid'):
        if(blockName in m_ElasticBlocksName):
            m_ElasticBlocksName.remove(blockName)  

def clearAllBlocksName():
    m_FluidBlocksName.clear()
    m_RigidBlocksName.clear()
    m_ElasticBlocksName.clear()

def printBlocksName():
    print('Fluids: ')
    for item in m_FluidBlocksName:
        print(item)
    print('RigidSolid: ')
    for item in m_RigidBlocksName:
        print(item)
    print('ElasticSolid: ')
    for item in m_ElasticBlocksName:
        print(item)