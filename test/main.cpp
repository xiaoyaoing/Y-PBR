
#include "Common/Render.hpp"
#include <iostream>

int main(int argc, const char* argv[]) {
    if (argc > 1) {
        for (int i = 1; i < argc; i++)
            Render::renderScene(argv[i]);
    } else {
        LOGE("No scene file provided please provide the folder containing the scene.json file as an argument.");
        Render::renderScene();
    }
}
