#include "src/encode4d.c"

#include "constants/4dDouble.h"
#include "ompExecBase.c"

int main()
{
  const struct CMUnitTest tests[] = {
    #include "testcases/omp.c"
  };

  return cmocka_run_group_tests(tests, setupRandomData, teardownRandomData);
}
