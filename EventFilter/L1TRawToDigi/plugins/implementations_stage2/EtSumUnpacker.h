#ifndef L1T_PACKER_STAGE2_ETSUMUNPACKER_H
#define L1T_PACKER_STAGE2_ETSUMUNPACKER_H

#include "EventFilter/L1TRawToDigi/interface/Unpacker.h"

namespace l1t {
   namespace stage2 {
      class EtSumUnpacker : public Unpacker {
         public:
            EtSumUnpacker();
            ~EtSumUnpacker() override {};

            bool unpack(const Block& block, UnpackerCollections *coll) override;

            inline void setEtSumCopy(const unsigned int copy) { EtSumCopy_ = copy; };

         private:
            EtSumBxCollection* res_;
            unsigned int EtSumCopy_;
      };
   }
}

#endif
