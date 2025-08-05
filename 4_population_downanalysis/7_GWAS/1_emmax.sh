#read -p "Input your plink_bedfile: " file_raw
#read -p "Input your plink_phefile: " phe_file
file=$1
phefile=$2
awk '{print $2}' $phefile > sample_phe.txt
awk '{print $2}' ${file}".fam" > sample_gen.txt

cat sample_phe.txt | xargs -n 1 -i grep -x {} sample_gen.txt > sample_gen_phe.txt
awk '{print 0, $1}' sample_gen_phe.txt > gen_phe.sample
plink --bfile $file --keep gen_phe.sample --make-bed --out ${file}.emmax --chr-set 27
rm sample*
plink --bfile ${file}.emmax --pca 10 --out ${file}.emmax --chr-set 27
cat ${file}.emmax.eigenvec | awk '{for(i=3; i<=NF; i++) printf $i " "; print ""}' > ${file}.emmax.gemma.covar 
awk '{print $1, $2, 1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ${file}.emmax.eigenvec > ${file}.emmax.emmax.covar
plink --bfile ${file}.emmax --recode12 --output-missing-genotype 0 --transpose --out ${file}.emmax --chr-set 27
awk '{print $2}' $file".fam" | xargs -n 1 -i grep -w '{}' $phefile > ${file}.emmax.phe
gemma -bfile ${file}.emmax -gk 2 -o ${file}.emmax -p ${file}.emmax.phe
mv output/* ./
for i in {3..14}; do echo "nohup gemma -bfile ${file}.emmax -n $i -k ${file}.emmax.sXX.txt -lmm 1 -o ${file}.emmax.${i}.covkin-TopTenPC -c ${file}.emmax.gemma.covar -p ${file}.emmax.phe &" | bash; done
for i in {3..14}; do echo "nohup gemma -bfile ${file}.emmax -n $i -k ${file}.emmax.sXX.txt -lmm 1 -o ${file}.emmax.${i}.covkin -p ${file}.emmax.phe &" | bash; done
emmax-kin-intel64 ${file}.emmax -v -d 10
emmax-kin-intel64 ${file}.emmax -v -s -d 10
for i in {3..14}
do
	awk -v i=$i '{print $1, $2, $i}' ${file}.emmax.phe > ${file}.emmax.tpm${i}.phe
        emmax-intel64 -v -d 10 -t ${file}.emmax -p ${file}.emmax.tpm${i}.phe -k ${file}.emmax.aIBS.kinf -o ${file}.emmax.${i}.assoc
        paste $file".emmax.bim"  ${file}.emmax.${i}.assoc.ps | awk '{print $1, $2, $4, $10, $7}' > ${file}.emmax.${i}.assoc.ps1
	mv ${file}.emmax.${i}.assoc.ps1 ${file}.emmax.${i}.assoc.ps
	emmax-intel64 -v -d 10 -t ${file}.emmax -p ${file}.emmax.tpm${i}.phe -k ${file}.emmax.aIBS.kinf -c ${file}.emmax.emmax.covar -o ${file}.emmax.${i}.covtopTenPC.assoc
	paste $file".emmax.bim"  ${file}.emmax.${i}.covtopTenPC.assoc.ps | awk '{print $1, $2, $4, $10, $7}' > ${file}.emmax.${i}.covtopTenPC.assoc.ps1
	mv ${file}.emmax.${i}.covtopTenPC.assoc.ps1 ${file}.emmax.${i}.covtopTenPC.assoc.ps 
	rm ${file}.emmax.tpm${i}.phe
done
