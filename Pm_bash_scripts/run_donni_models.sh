cd "/mnt/c/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun"

for i in bottlegrowth_1d growth snm_1d three_epoch three_epoch_interbreeding two_epoch

do donni infer --input_fs dadi_1pop_wsaf.fs --model ${i} --export_dadi_cli 1D_demographics/${i} --output_prefix 1D_demographics/${i}

done

#for i in IM IM_mscore IM_pre IM_pre_mscore bottlegrowth_2d bottlegrowth_split bottlegrowth_split_mig snm_2d split_asym_mig split_delay_mig split_mig split_mig_mscore no_mig sym_mig asym_mig no_mig_size sym_mig_size asym_mig_size sec_contact_sym_mig sec_contact_asym_mig anc_sym_mig anc_asym_mig sec_contact_sym_mig_size sec_contact_asym_mig_size anc_sym_mig_size anc_asym_mig_size sym_mig_twoepoch asym_mig_twoepoch sec_contact_sym_mig_three_epoch sec_contact_asym_mig_three_epoch sec_contact_sym_mig_size_three_epoch sec_contact_asym_mig_size_three_epoch

#do donni infer --input_fs dadi_admixture_no_af_filter.fs --model ${i} --export_dadi_cli 2D_demographics/${i} --output_prefix 2D_demographics/${i}
#done


#for i in split_nomig split_symmig_all split_symmig_adjacent ancmig_adj_3 ancmig_adj_2 ancmig_adj_1 sim_split_no_mig sim_split_no_mig_size sim_split_sym_mig_all sim_split_sym_mig_adjacent split_nomig_size ancmig_2_size split_sym_mig_adjacent_var1 split_sym_mig_adjacent_var2 split_uni_mig_adjacent_var1 split_uni_mig_adjacent_var2 sim_split_sym_mig_adjacent_var sim_split_uni_mig_adjacent_var admix_origin_no_mig admix_origin_sym_mig_adj admix_origin_uni_mig_adj

#do donni infer --input_fs dadi_DAPC_no_af_filter.fs --model ${i} --export_dadi_cli 3D_demographics/${i} --output_prefix 3D_demographics/${i}
#done
