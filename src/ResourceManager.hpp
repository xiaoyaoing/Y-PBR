#pragma once

#include "memory"

class Scene;
class Integrator;
class Camera;
class Image;

class ResourceManager {
public:
    static std::shared_ptr<Scene> getScene() {
        if (!_scene.get()) {
            throw("Scene not initialized!");
        }
        return _scene;
    }

    static void setScene(std::shared_ptr<Scene> scene) {
        _scene = scene;
    }

private:
    static std::shared_ptr<Scene>      _scene;
    static std::shared_ptr<Integrator> _integrator;
    static std::shared_ptr<Camera>     _camera;
    static std::shared_ptr<Image>      _image;
};