# ![](logo.png) SPHinXsys Refactoring 

## List of smells

1. it seems that there is a bug for reading stl file (the file name issue?) when using debug mode in windows 
2. Introduce boundary condition or constraints which go together with particle dynamics algorithms. In this way, user will not forget the boundary condition or constraints during time stepping? 
3. naming issues: used base key (like base_particles_) in name when it really means base class objects, if virtual function overriding is expected, do not use it (like use particles_)
4. closest point on surface for 3d geometries with binary operations
5. using more meaningful names for class members, so that I do not need extra comments
6. xml memory leaking (tried second time. Failed. Erasing XML elements will not decrease the memory usage.)
7. user may forget impose constraint when apply an algorithm modifying velocity field
8. temporary and one time output, such as for level set, may open a temporary folder without related to system.  

## General principles on reviewing the code

SPHinXsys should be reviewed often to increase its quality continuously. The target is to achieve a good layout which arranges the corresponding classes, variables, functions, algorithms, methods and physics in reasonable way. Generally, there are several criteria, i.e. formality, locality, additivity, simplicity and consistency, to qualify the layout of a code.  

### Formality

As a C++ code, we would like to follow [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

### Locality

First, locality requires that the data to operated should be as close to each other as possible. Second, it requires that conceptionally related variables and methods should be located close to each other.

### Additivity

When a new feature is added, the original functionality should be maintained in general. This may requires modification of the original code. However, such modification should not change the interface by which the original functions are called, except that the interface will also be updated.

### Simplicity

Less lines of code in a user case indicates easier application of SPHinXsys. A design will be more attractive if it can also avoid mistakes when a user case is set up. There are several ways to increase the simplicity. One is to define several relevant parts once for all.

### Consistency

Consistency is very important. There are many aspects of consistency. The name of variable, function or class should be meaningful and consistent with its functionality. The comment should be consistent with the actual usage of the corresponding entry.

## Narrative coding

One idea is writing a code in narrative style, just like to write a story. Since the code is quite complex, so it is like a complex story (may be a computer game?) with many stroy lines.

## Error-driving optimization

After an error in an application or the library is reported, when we try to correct the error, we need to think what is possible reason that the error occurs, and try to optimize the code for avoiding in future similar errors.

## User-experience emphasis

When we need choose between small performance gain and code tolerance for less strict coding style and better user experience, we more likely choose for the latter.