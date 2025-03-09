tag=$1


ANNOVARDIR=../opt/annovar
OUTPUTDIR=../data/processed/annovar
TARGETDIR=../data/processed/pgs


cat ${TARGETDIR}/${tag}_pgs_weights_hg38.txt \
	| tail -n +2 \
	| awk -v OFS='\t' '{print $1,$2,$2,$5,$4}' \
	> annovar_input.tmp

perl $ANNOVARDIR/table_annovar.pl \
	annovar_input.tmp \
	${ANNOVARDIR}/humandb/ \
	-buildver hg38 \
        -out ${OUTPUTDIR}/${tag} \
	-protocol refGene \
	-operation g \
	-remove -polish -nastring .

rm annovar_input.tmp
