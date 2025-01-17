# Y-PBR
Physically Ray Tracer 

Integrators:
- Path Tracing
- Bidirectional Path Tracing
- SPPM
- Volume Path Tracing

Materials:
- Lambertian
- Mirror
- Dielectric
- Conductor
- Plastic
- Hair

Volume:
- Homogeneous

Gallary
![kitchen.png](Gallery%2Fkitchen.png)
![ajar.png](Gallery%2Fveach-ajar.png)
![bidir.png](Gallery%2Fbidir.png)![veach-bidir12-4.png](Gallery%2Fsppm.png)
![bidir-water.png](Gallery%2Fbidir-water.png)
![classroom.png](Gallery%2Fclassroom.png)
![classroom1.png](Gallery%2Fclassroom1.png)
![star.png](Gallery%2Fstar.png)
![star1.png](Gallery%2Fstar1.png)

Reference

Scene Mainly From https://benedikt-bitterli.me/resources/

Algorithm Reference From https://www.pbrt.org/  


Tungsten Render  https://github.com/tunabrain/tungsten  

Nori Render https://wjakob.github.io/nori/  

Ray Tracing in One Weekend https://raytracing.github.io/books/RayTracingInOneWeekend.html  


Build
Run the following command in the root directory of the project
```shell
mkdir build
cd build
cmake -G "Visual Studio 17 2022" ..
```
Then open the solution file in the build directory and build the project.

Usage:
```shell
Y-PBR.exe scene.json
```
[scene.json](scenes/classroom/scene.json)
Scenes can be found in the Scenes directory.
