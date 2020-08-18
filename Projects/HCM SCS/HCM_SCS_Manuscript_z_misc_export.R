


# Export the count tables only for purposes of comparison
hcm_scs = groupedSCS$patientAllMod
save(list = 'hcm_scs', file = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/_sessions/export_SCS_only.Rdata')
rm('hcm_scs')


# Export PKP2 vulcano
PKP2_vol_df = get_volcano_df(as.matrix(groupedSCS$patientAllMod@ndata), 'PKP2', calc_qvals = F)
plot_volcano(PKP2_vol_df, NRLABELED = 50)+xlim(c(-.25,.25))+ggtitle('Genes correlated with PKP2')
