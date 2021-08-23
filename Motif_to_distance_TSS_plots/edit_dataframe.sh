for i in $@
do
Rscript edit_dataframe_for_plot.R $i $i.bed
done
