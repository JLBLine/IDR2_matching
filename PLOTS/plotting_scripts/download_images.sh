# python download_images.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt -c ../../multifreq/gleammulti_reject_sources.txt
# python download_images.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt -c ../../multifreq/gleammulti_eyeball_sources.txt
# python download_failed.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt

python download_images.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt -c ../../IDR2_eyeball_sources.txt
python download_images.py -i ../../deep/puma_gleamdeep-v-m-s-n-a-eyeball.txt -c ../../IDR2_eyeball_sources.txt
python download_images.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt -c ../../IDR2_reject_sources.txt
python download_images.py -i ../../deep/puma_gleamdeep-v-m-s-n-a-eyeball.txt -c ../../IDR2_reject_sources.txt
# python download_failed.py -i ../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt
python download_failed.py -i ../../deep/puma_gleamdeep-v-m-s-n-a-eyeball.txt