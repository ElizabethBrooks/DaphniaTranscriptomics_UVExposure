# script to subset mapping stats by sample group

# subset olympic

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/alignment_hisat2_jobOutput.o* \
| grep -B1 -e "MUV10_" -e "MUV27_" -e "MUV28_" -e "MUV17_" -e "MUV20_" -e "MUV3_" -e "MUV31_" -e "MUV32_" -e "MUV33_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/alignment_hisat2_EGAPx_overall_olympic.txt

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/alignment_hisat2_test_jobOutput.o* \
| grep -B1 -e "MUV10_" -e "MUV27_" -e "MUV28_" -e "MUV17_" -e "MUV20_" -e "MUV3_" -e "MUV31_" -e "MUV32_" -e "MUV33_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/alignment_hisat2_EGAPx_test_overall_olympic.txt

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/alignment_hisat2_test_jobOutput.o* \
| grep -B1 -e "MUV10_" -e "MUV27_" -e "MUV28_" -e "MUV17_" -e "MUV20_" -e "MUV3_" -e "MUV31_" -e "MUV32_" -e "MUV33_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/alignment_hisat2_test_overall_olympic.txt

# subset sierra

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/alignment_hisat2_jobOutput.o* \
| grep -B1 -e "MUV11_" -e "MUV30_" -e "MUV6_" -e "MUV13_" -e "MUV16_" -e "MUV1_" -e "MUV24_" -e "MUV25_" -e "MUV7_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/alignment_hisat2_EGAPx_overall_sierra.txt

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/alignment_hisat2_test_jobOutput.o* \
| grep -B1 -e "MUV11_" -e "MUV30_" -e "MUV6_" -e "MUV13_" -e "MUV16_" -e "MUV1_" -e "MUV24_" -e "MUV25_" -e "MUV7_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/alignment_hisat2_EGAPx_test_overall_sierra.txt

cat /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/alignment_hisat2_test_jobOutput.o* \
| grep -B1 -e "MUV11_" -e "MUV30_" -e "MUV6_" -e "MUV13_" -e "MUV16_" -e "MUV1_" -e "MUV24_" -e "MUV25_" -e "MUV7_" \
| grep "overall" | sed "s/% /,/g" > /Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/alignment_hisat2_test_overall_sierra.txt
