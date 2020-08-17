# ![](logo.png) Coding SPHinXsys

## To do list

1. revise class members by checking Doxgen document
2. Niklas Aufenanger's case goes strange (bubbles at the wall boundary) when resolution is very high.
3. body update state
4. closest point on surface for 3d geometries with binary operations
5. direction graph for linear shapes in a SPH body
6. xml memory leaking (tried second time. Failed. Erasing XML elements will not decrease the memory usage.)
7. user may forget impose constraint when apply an algorithm modifying velocity field
8. FSI when consider more complex wall models for high Reynolds number flows
9. using more meaningful names for class members, so that I do not need extra comments
10. resolving the compiling warnings for the order of class member initialization
11. lower- and upper bounds can be represented by a pair

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