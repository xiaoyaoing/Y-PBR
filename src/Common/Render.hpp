#pragma once

#include "Integrator/PathIntegrator.hpp"
#include "work-queue.hpp"
#include "IO/FileUtils.hpp"

struct Bucket {
    Bucket() : min(0), max(0) {}
    Bucket(const ivec2& min, const ivec2& max) : min(min), max(max) {}

    ivec2 min;
    ivec2 max;
};

class Render {
public:
    static void renderScene(const std::string path = FileUtils::WorkingDir);

    Render(const Json& json);
    void Go();

public:
    std::shared_ptr<Camera>     camera;
    std::shared_ptr<Scene>      scene;
    std::unique_ptr<Integrator> integrator;
    std::shared_ptr<Sampler>    sampler;
};

class TimeCounter {
public:
    void printTimeCount() {
        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        auto                                               s   = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        std::cout << std::endl
                  << event << "done. Take " << s << "s";
    }

    template<typename type = std::chrono::milliseconds>
    static void timeClock(long long time) {
        TimeCounter counter;
        while (counter.getTimeCount<type>() < time)
            ;
    }
    template<typename type = std::chrono::milliseconds>
    long long getTimeCount() {
        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        return std::chrono::duration_cast < type >> (end - start).count();
    }
    TimeCounter(std::string event = "") : event(std::move(event)) {
        start = std::chrono::steady_clock::now();
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::string                                        event;
};