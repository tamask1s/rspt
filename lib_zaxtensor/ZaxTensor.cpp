/***************************************************************************
* Copyright 2020 Tamas Levente Kis                                         *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/

#include <typeinfo>
#include <string>
#include <string.h>
#include <map>
#include <vector>
#include <numeric>
#include <stdio.h>
#include <cstdlib>

using namespace std;

#include "ZaxJsonParser.h"
#include "ZaxTensor.h"

vector<int> get_dimensions(const char* a_json)
{
    vector<int> res;
    const char* a_str_to_find = a_json;
    int nrdimensions = 0;
    while (*a_str_to_find)
        if (*a_str_to_find++ == '[')
            ++nrdimensions;
        else if (*a_str_to_find == ']')
            break;

    res.resize(nrdimensions);
    a_str_to_find = a_json;
    int deepness = -1;
    while (*a_str_to_find)
        if (*a_str_to_find++ == '[')
            ++deepness;
        else if (*a_str_to_find == ',')
            ++res[deepness];
        else if (*a_str_to_find == ']')
            ++res[deepness--];
    for (unsigned int i = res.size() - 1; i > 0; --i)
        res[i] /= res[i - 1];

    return res;
}
