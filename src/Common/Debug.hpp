//2022/8/26
#pragma once

struct DebugConfig{
    const static  bool OnlyDirectLighting;
    const static  bool OnlyIndirectLighting;
    const static  bool OnlyOneThread;
    const static  bool OnlyShowNormal;
    const static  bool DebugMode ;


    const static  bool sampleInRange = false;
    const static  int  sampleMinx;
    const static  int  sampleMiny;
    const static  int  sampleMaxx;
    const static  int  sampleMaxy;
};
