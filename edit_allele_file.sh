#!/bin/bash
mv trimmed.tfa emm_tsdna.fas
sed -i -e 's/EMM/emm/g' emm_tsdna.fas
sed -i -e 's/ emm/.sds emm/g' emm_tsdna.fas
sed -i -e 's/  ASSEMBLE/.sds  ASSEMBLE/g' emm_tsdna.fas
sed -i -e 's/ st/.sds st/g' emm_tsdna.fas
