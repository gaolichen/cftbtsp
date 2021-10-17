#include <boost/test/unit_test.hpp>
#include "test.h"

std::vector<CftConfig> CreateCfgConfigTestData(int size)
{
    std::vector<CftConfig> ret;
    for (int i = 0; i < size; i++) {
        CftConfig config(randomint(2, 7));
        int maxSpin = randomint(3, 5);
        for (int j = 0; j <= maxSpin; j++) {
            config.OperatorNumbers.push_back(randomint(3, 5));
        }

        ret.push_back(config);
    }
    
    return ret;
}

