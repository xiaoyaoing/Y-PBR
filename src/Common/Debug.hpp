//2022/8/26
#pragma once

struct DebugConfig{
    const static  bool OnlyDirectLighting= false;
    const static  bool OnlyIndirectLighting= false;
    const static  bool DebugMode ;
    const static  bool OnlyOneThread= false;
    const static  bool OnlyShowNormal  =false;

    const static  bool sampleInRange = false;
    const static  int  sampleMinx;
    const static  int  sampleMiny;
    const static  int  sampleMaxx;
    const static  int  sampleMaxy;
};
