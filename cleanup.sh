#!/bin/bash
ls | grep -xv "produce_edep-paramreco_muoncontained.sh\|ND_CAFMaker_job_uptoedep_nogapdset\|cleanup.sh" | xargs rm -r
