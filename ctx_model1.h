/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
* Author(s): Marcus D. R. Klarqvist
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*/

#include "pbwt.h"

namespace djinn {

class djinn_ctx_model {
public:
    int Compress(uint8_t* data);

    std::shared_ptr<GeneralModel> mref;
    std::shared_ptr<GeneralModel> mlog_rle, mlog_rle_o1;
    std::shared_ptr<GeneralModel> mrle, mrle_o1;
    std::shared_ptr<GeneralModel> mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4;
};

}