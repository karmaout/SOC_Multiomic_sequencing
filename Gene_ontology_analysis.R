library(clusterProfiler)
library(org.Rn.eg.db)
library(DOSE)
library(ggplot2)
library(openxlsx)
library(stringr)
library(readxl)

# Function to perform GO analysis and return results as a data frame
perform_go_analysis <- function(gene_list, group_name) {
  # Convert gene symbols to Entrez IDs
  gene_df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
  
  # Identify unmapped genes
  unmapped_genes <- setdiff(gene_list, gene_df$SYMBOL)
  cat("Unmapped genes for", group_name, ":", unmapped_genes, "\n")
  
  # Use only the mapped genes for GO analysis
  mapped_genes <- gene_df$ENTREZID
  
  # Perform GO enrichment analysis
  ego <- enrichGO(gene = mapped_genes,
                  OrgDb = org.Rn.eg.db,
                  ont = "BP",  # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.2,  # More lenient cutoff
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  # Convert results to data frame
  ego_df <- as.data.frame(ego)
  
  # Add group name to the results if there are any results
  if (nrow(ego_df) > 0) {
    ego_df$Group <- group_name
  } else {
    ego_df <- data.frame()
  }
  
  return(ego_df)
}

# Define gene lists for each group
all_DEG_lists <- list(
  VGLUT1_UE = c(
    "Pou3f1", "Haao", "Fut7", "Depdc7", "Top2a", "Ccdc42", "Sp6", 
    "Calm-ps2", "Atp13a5", "Ubxn2a.1", "H1f10", "Pvr", "Ttc39d", 
    "Prss30", "Amtn", "Haus4", "Sik1", "Sh2b2", "Ervfrd-1", "Sh3bp1", 
    "Ceacam9", "Oaz3", "Dlec1", "Rec8", "Vstm4", "Adamts12", "Fos"
  ),
  VGLUT1_PE = c(
    "Fgf21", "Meiosin", "Olr331", "Nlrp5", "Vom1r100.1", "Ankrd30a", 
    "Slc12a1", "AC131483.1", "Twist1", "Otogl", "Crb1", "Dpm3.1", 
    "Ecscr", "Olr1668", "Nr4a2", "Arc", "Trib1", "Egr1", "Pipox", 
    "Igsf10", "Vom1r46", "Bnc1", "Vom2r23", "Chrna10", "AC096024.2", 
    "Cd1d1", "Fcrl2", "Ccdc141", "Rpl32-ps3", "Bpifb4", "Nrarp", 
    "Ttc30a1", "Tmem200b", "Dnaaf6", "Inhbc", "Snx33", "Ggt6", 
    "Cd300lg", "Alox15", "Rab7b", "Slc26a9", "Tmem253", "Myh6", 
    "Nid2", "Eci3", "RGD1565323", "Olig2", "NEWGENE-2724", "Cyyr1", 
    "Hsf2bp", "Nr4a1", "Cavin3", "Homer1", "Egr4", "Pef1", "Ednrb", 
    "Vom1r47", "Pard6a", "Fos"
  ),
  VGLUT1_UR = c(
    "Tiparp", "Egr2", "Nr4a2", "Nr4a1", "Egr4", "Egr1", "Adrb1", 
    "Abra", "Adcyap1", "Bach1", "Bdnf", "Btg1", "Rasl11a", "Commd2", 
    "Homer1", "Dusp6", "Thrb", "Syngr4", "Serpinh1", "Cox7a1", "Itgam", 
    "RGD1562492", "Lipm", "Hamp", "As3mt", "Samd3", "Cep55", "Ccdc172", 
    "Gng5", "Pdlim5", "Clec2e", "Gprc5a", "Trib3", "Spata25", "Defb27", 
    "Olr546", "Crispld1", "Tnfrsf14", "Cyp2j16", "Gpr3", "St18", "Esrp1", 
    "Atp1b4", "Shisa8", "RGD1561102", "RGD1566368", "Tm4sf20", "Arid5a", 
    "Il18rap", "Wdr64", "Ptprv", "Hist1h2bo", "Adamts1", "Sik1", "Ephb4", 
    "Ptgs2", "Fos"
  ),
  VGLUT1_PR = c(
    "Arc", "Bdnf", "Zfp773-ps1", "Pde7b", "Scg2", "Elmo1", "Cd2", 
    "Ptx3", "Gprc5a", "Bmp2", "Bdkrb2", "Klhl30", "Fam178b", "Abhd15", 
    "Kcnmb1", "Stc2", "Irgm", "Zp3r", "NEWGENE-1310011", "Gpa33", 
    "Rab5al1", "Atp13a4", "Ttc29", "RT1-DOb", "Ntrk2", "Tiparp", 
    "Egr1", "Pde10a", "Nptx2", "Trib1", "Egr4", "Plcxd2", "Nr4a1", 
    "Lingo1", "Poted", "Dusp6", "Sorbs2", "Kcnk3", "Acot3", "Tmc8", 
    "Arl9", "Ncam1", "Ddr2", "Snap25", "Asb12", "Htra4", "Slc35e1", 
    "Mapk4", "Maml3", "RGD1563072", "Lxn", "Ackr1", "Fos"
  ),
  Sst_PE = c(
    "Nfkbiz", "Vom2r38", "Hapln3", "Vamp8", "Cfap157", "Tal2", 
    "Col9a2", "Radx", "Arid5a", "Hes7", "Olr1410", "Lax1", "Lgals3", 
    "Cyp4v3", "Sncg", "Nptx2", "RGD1564941", "Mtmr11", "Spsb3", 
    "Tent5a", "Vgf", "Spred3", "Nr4a1", "Tln1", "Tmie", "Aqp4", 
    "Slc5a4", "Gadd45g", "Dmgdh", "Kcns3", "Fcgbpl1", "Ppp1r3c", 
    "Fgl2", "Sinhcaf", "Fbln2", "Olfml2a", "Shisal2a", "Col16a1", 
    "S100g", "Abcb5", "Acp5", "RGD1560917", "NEWGENE-620180", 
    "RGD1561662", "Tcap", "Ncf2", "Tex35", "Dpm3.1", "Mcm10", "Casr", 
    "Itih4", "Tnfaip1", "Snrpa", "Uap1l1", "Gpr6", "Tspyl5", "Fos"
  ),
  Sst_UR = c(
    "Prok2", "Esrra", "Ankrd1", "Irag1", "Vcam1", "Nudt17", 
    "Hapln2", "Pon1", "Il23r", "Sinhcaf", "Il5ra", "Ccdc183", 
    "Wfdc10a", "Spaca1", "Magee2", "Vsx2", "Cyp4f6", "Trib1", 
    "Foxh1", "Cxcr5", "Kcnj1", "Snorc", "Igfbp4", "Abhd15", 
    "Fam163a", "Pla2g4a", "Tmem273", "RGD1561777", "Clps", 
    "RT1-CE16", "Clec4m", "Kcnj16", "Vgf", "Shld3", "Slc25a18", 
    "Cdhr1", "Cnppd1", "RGD1563888", "Efcab3", "Tomm40l", 
    "Hopx", "Adra1b", "Myh7", "Cbfb", "Wdcp", "Gpr160", 
    "Pdlim5", "Inhba", "Eya4", "Heyl", "Aqp9", "Prrx1", 
    "Mesp2", "Gldc", "Hapln3", "Ptprcap", "Afap1l2", 
    "Mpeg1", "Dnase2b", "Vtcn1", "Ankdd1b", "Wnt2", "Serinc4", 
    "Upp2", "Helz2", "Rex2", "Col16a1", "Calr4", "Xkrx", 
    "Rtl1", "Spatc1", "Ooep", "Mstn", "Twist2", "Aox4", 
    "Ctla4", "Gp1ba", "Mael", "Spata46", "Scel", "Runx1", 
    "Apela", "Plvap", "Pcdhb7", "Slc15a3", "Stom", 
    "Ccdc105", "Naprt", "Bace2", "Scamp5", "Fgf2", 
    "Mtmr11", "Tomm5", "Yod1", "Rgs20", "Znfx1", 
    "Ipcef1", "Rwdd4", "Fos"
  ),
  Sst_PR = c(
    "Plg", "Il33", "Vtcn1", "Ccn1", "Itga1", "Dmrta2", "Mybpc1", 
    "Als2cl", "Stbd1", "H1f1", "RT1-Db1", "Trib1", "Amotl2", 
    "Afap1l2", "Aldh1a2", "NEWGENE-620180", "Serpinf2", "Mettl11b", 
    "Prom1", "Olig2", "AC110981.1", "Ppic", "Gpr17", "Ccdc178", 
    "Cplx4", "Chtf18", "Ednrb", "C2cd4d", "Nr4a1", "Mfsd10", 
    "Trim5", "Olr319", "Rab13", "Tril", "Knl1", "Tie1", 
    "Sap18", "Lrrc63", "Nr2e1", "RGD1561426", "Etv5", "Kif4b", 
    "Tp53inp2", "Phlda1", "Trip13", "Lair1", "Il12a", "Lxn", 
    "Gsg1", "Col27a1", "Rp1", "Tbxa2r", "Olr1307", "C2cd4b", 
    "Scn10a", "Efcab3", "Rps18l1", "Lmod1", "Fmod", "Nkapl", 
    "Slc15a2", "Adam21", "Cdkn1a", "Matn4", "Nr4a2", "Pla2r1", 
    "Lingo3", "Pdgfd", "Slc35g3", "Fam89a", "Cspg4", "Tmem200a", 
    "Sbsn", "Il4r", "Lurap1l", "Serpind1", "Galk1", "Hus1", 
    "Arl10", "Abhd13", "Styk1", "Cmpk2", "Abcb1b.1", "Trpv6", 
    "Tmem269", "Samd11", "Zfp36l2", "Abcb5", "Gpr18", "Upb1", 
    "Tnip1", "Slc5a2", "Cdc42ep4", "Bambi", "Mboat4", "Csdc2", 
    "Rbm48", "R3hdm4", "MGC116202", "Vps18", "C1ql1", "Nfkbia", 
    "RGD1562402", "Grem1", "Bub1", "Slco5a1", "Tmem45b", "Amigo3", 
    "B3gnt7", "Catsper3", "Tpm4", "Cyba", "Ostf1", "Tinagl1", "Fos"
  ),
  MSN_UE = c(
    "Papln", "Bub1", "Tns1", "Odad1", "Abcg5", "Olr93", "AC118496.1", 
    "Fbln2", "Btg2", "Gdpd3", "Il27", "Cd2", "Pdlim5", "Chrm2", 
    "Slc41a3", "Rpl32-ps3", "Pla2r1", "Pamr1", "4930444P10Rik", 
    "Megf6", "Tnfrsf9", "Oxct2a", "Gem", "Tbx22", "Cyp4f6", "Igfbp6", 
    "Tg", "Epor", "Snrpel1", "Impdh2", "Izumo1r", "Znrd1as", 
    "RGD1566007", "Dnajb3", "Arhgap28", "Glis2", "Ciita", "Dtl", 
    "Tug1", "Gck", "RGD1308117", "Tdh", "Slc51a", "Cyp4f18", 
    "Csf1r", "Slc6a7", "Palm3", "Ucp1", "Aif1", "RGD1561777", "Cdkn1a", "Fos"
  ),
  MSN_PE = c(
    "Syngr4", "Siglec5", "Cd37", "Vwce", "Fes", "Lilrb3a", "Tmigd3", 
    "Mex3a", "Enpep", "Clec2g", "Bcl2l11", "Spaca9", "Pamr1", 
    "Gna15", "Bmp5", "Ankrd23", "Card14", "Galr2", "Tk1", 
    "B4galnt2", "Slamf1", "RGD1310587", "Cenpf", "Rgs1", 
    "Fcgr2b", "Slc10a6", "Map3k8", "Ifi30", "Osgin1", "Fgd2", 
    "Slc19a1", "Pus7l", "Myoz2", "Arhgdib", "Ptprc", "Ms4a3", 
    "Trim5", "Tkfc", "Aknad1", "Cd101", "Dbf4", "Ptpn6", 
    "Xirp2", "Tmem269", "Shisa8", "Fli1", "Sp140", "Nuf2", 
    "Sh3bp2", "Olr1653", "Chrd", "Itgb5", "Il13ra1.1", "Pros1", 
    "Stard4", "Rnaseh2a.1", "RT1-Ba", "Icoslg", "Clec4a3", 
    "Col11a2", "Dnah8", "Pcsk1", "Nrg4", "Rabl2", "Oasl2", 
    "Trim59", "Klra2", "Spi1", "Asb12", "Gja6", "RGD1565071", 
    "Fut4", "Pkhd1", "Ltc4s", "Pik3r6", "Mmgt2", "Hmmr", 
    "Hmcn1", "Ska3", "Syk", "Runx1", "Myo18b", "Pgr15l", 
    "Pdp2", "Slc7a8", "Fos"
  ),
  MSN_UR = c(
    "Csf1r", "Endog", "Ddx28", "Slco2b1", "Csrnp1", "Acer1", 
    "Igfbp4", "Synb", "Plekhm2", "Skiv2l", "Cry1", "Afap1", 
    "Arv1", "Slc45a2", "Adra2b", "Nfatc2", "Yipf7", "Mtres1", "Fos"
  ),
  MSN_PR = c(
    "Pipox", "Adam5", "C5ar2", "Gnmt", "Nek2l1", "Thbs2", "Trim34", 
    "Pgghg", "Sipa1", "Cyp2c23", "Vom2r16", "Mcub", "Emx1", 
    "Dppa3", "Tlr4", "Samd11", "RGD1560314", "Vom2r57", "Casp1", 
    "Sln", "Ccr3", "Olr1516", "Kif19", "P2rx1", "Abi3", 
    "RGD1564463", "Atf3", "Fam114a1", "Fgfr3", "Rnf212b", 
    "Hist1h2ail1", "Primpol", "Psma8", "Sall3", "Ier2", "Nr2e1", 
    "Sapcd1", "Sgo2", "Cc2d1a", "Arhgap27", "Sart1", "Sfmbt2", 
    "Arc", "RT1-T24-1", "Selenop", "Tspan8", "Nde1", "Lacc1", "Adgrd1", "Fos"
  ),
  NPY_PE = c(
    "Slc16a12", "Igfbp2", "Bcl2l12", "MGC94891", "Cacng6", "Dnase2b", 
    "P2ry13", "Dbf4", "Thbd", "Scn7a", "Rbm38", "RGD1306233", 
    "Efcab11", "Six4", "Maff", "Hmga2", "Mei1", "Crtam", "Minar1", 
    "Il20rb", "Col5a2", "Slc11a1", "Rab37", "Usp43", "Elf3", 
    "Cacna1s", "Avpr1b", "Ikbke", "Gck", "Gch1", "Hes1", "Hmox1", 
    "Cd248", "Ltbp1", "Baiap2l2", "Ptchd4", "Etv4", "Qpctl", 
    "Trib1", "Serpind1", "Ccdc159", "Man1a1", "Bri3bp", "Parp11", 
    "Ntn4", "Slc7a7", "Trip6", "Irak4", "Phex", "Vmac", "Taf1c", 
    "Syt10", "Sim2", "Sox6", "Lrrc27", "Arnt2", "Pgm5", "Moxd1", 
    "RGD1559962", "Myof", "Lrrc70", "Amy1a", "Prrg4", "Tgm7", 
    "E2f2", "Ube2u", "Heph", "Angpt1", "Fbxo32", "Stac", "Trip10", 
    "Obscn", "Sstr2", "Grm6", "Plcd3", "Arhgap27", "C1qtnf7", 
    "S1pr3", "Fgf20", "Slc8b1", "Pole", "Tcf7l1", "Col5a3", "Fos"
  ),
  NPY_UR = c(
    "Susd2", "Sik1", "Mamdc2", "Man1a1", "Vsig8", "Ephb3", "Asphd1", 
    "Btg2", "Cbx8", "Tbc1d10c", "Klk6", "Ano9", "Fam53b", "Fpr1", 
    "Ahrr", "Tbx10", "Slc35g1", "Bglap", "Dennd2d", "Mab21l1", 
    "S100a10", "Atoh8", "Herc6", "B4galnt3", "Abo3", "Ripor3", 
    "Mybl1", "Cavin4", "Timp1", "Rps10l1", "Zfp36l1", "Foxh1", 
    "Gpr62", "Bard1", "Sema4c", "Hes7", "Acap1", "Efcab3", 
    "AC130970.1", "Cby3", "Tnn", "Sh3tc1", "Hist2h3c2.3", "Card19", 
    "Hes1", "Adprhl1", "Megf10", "Cd74", "Spry4", "Syce1l", 
    "Nudt7", "Dnase2", "Naglt1", "Rhobtb1", "Pcolce", "Chek2", 
    "Tmem233", "Golga3", "Brcc3", "Bambi", "Fos"
  ),
  Npy_PR = c(
    "Nr4a1", "Slc41a3", "Cnksr3", "Arc", "Egr1", "Fblim1", "Chac1", 
    "Nr4a3", "Trib1", "Vgf", "Man1a1", "Etv5", "Itprip", "Igsf23", 
    "Plin1", "1700092M07Rik", "AC098622.1", "Tbx10", "Crh", "Gpr160", 
    "Spata16", "Rab13", "Them5", "Zfp862", "C5", "Lcn2", "Pdrg1.1", 
    "Gpr3", "RGD1564855", "Pgr15l", "Ism2", "Dio2", "Avpr1a", 
    "Iqcf1", "Gcnt3", "Paqr5", "Ccdc69", "Krt27", "Prr29", 
    "Tnfsf18", "Cysltr2", "Olr1658", "Kng2l1", "Cmtm2a", "Bean1", 
    "Klf1", "Nlrc5", "Olr1700", "Rxfp2", "Fgf1", "Stimate", 
    "Lrrtm3", "Mkx", "Kat2b", "Hs3st4", "Taf5", "Cacng4", 
    "Pard6a", "Letmd1", "Pomc", "Sgcz", "Egr2", "Gnb5", "Creb3l1", 
    "Pamr1", "Gpr83", "RGD1561777", "Adora2a", "Sptssb", 
    "Marcksl1", "Gramd1a", "Necab1.1", "Siah2", "Nebl", "Frat2", "Rasgrp2", "Fos"
  ),
  PV_UE = c(
    "RGD1561413", "Cnppd1", "Ninj2", "Bcas1", "Rmc1", "Agmo", 
    "Mobp", "Synm", "Nudt8", "Apip", "Arrdc1", "Ctdsp1", 
    "Dhrs11", "Cacng5", "Synpo2l", "Fos"
  ),
  PV_PE = c(
    "Wdr49", "Cpvl", "Emx1", "Ttll9", "Tp53inp1", "RGD1565410", 
    "Fbxw12", "C2cd4a", "Fbxl22", "Nicn1", "Jmjd4", "Nr1i3", 
    "Cxcl6", "Dmp1", "Pdgfrb", "Vhl", "Mrpl10", "Hyi", 
    "Senp3", "Mepce", "Fpr1", "Nfkbid", "Arhgef38", "Tex48", 
    "Dnph1", "Mslnl", "Nsd3", "Cibar2", "Tspyl1", "Chst9", "Fos"
  ),
  PV_UR = c(
    "Glra1", "Ftsj3", "RGD1562229", "Lrrc36", "Cdc25b", "Nkd2", 
    "Prrg2", "Omp", "Tspan32", "Relb", "Cdkn1c", "Rln1", 
    "Marveld1", "Calca", "Tnfsf10", "Lrrc71", "Ptger3", "S100a1", 
    "Arsj", "Ret", "Wnt2", "Svopl", "Dyrk4", "Vom1r90", 
    "Fbln2", "Prg2", "Rassf2", "Adamtsl2", "Akr7a3", "Col15a1", 
    "Cyb5rl", "Stpg4", "Smoc1", "Pdgfb", "Slc17a8", "Gm9918", 
    "Tamalin", "Ephx3", "Rasl12", "Clec3b", "Tnfsf9", "Mreg", 
    "Trip10", "Trpv1", "Cluap1", "Abca6", "Per1", "Meikin", 
    "Nuak2", "Etnk2", "Lhx4", "Rps29.1", "Jchain", "C1qtnf7", 
    "Mcpt4l1", "Dok3", "Hist2h3c2.4", "Fam217a", "Fgl1", "Plac8l1", 
    "Ccl17", "Olr1694", "Tysnd1", "RT1-DMb", "Ercc2", "Egr4", "Elk4", "Fos"
  ),
  PV_PR = c(
    "Hdc", "Bbc3", "Ccdc96", "Colgalt1", "Steep1", "Med7", "Tshz3", 
    "Calhm2", "Mmp21", "Ctsc", "RGD1563307", "Hey2", "Ms4a2", 
    "Sycp1", "Tshb", "Flnc", "Chmp4bl1", "Kynu", "AC127963.2", 
    "Ttc30a1", "Tldc2", "Ctrc", "Serpina11", "Cyth4", "Cthrc1", 
    "Sla", "Notch3", "Pate1", "AC128962.3", "Tram2", "Pimreg", 
    "Pkmyt1", "P2rx1", "Atp5mk.1", "Gast", "Lad1", "Pappa2", 
    "Pla2g3", "Fam3d", "Hist1h2an", "Itih4", "Lox", "Myo7b", 
    "Pde6a", "Il27ra", "Crispld2", "Ces1c", "Spc24", "Ccnf", "Marchf8", "Fos"
  ),
  Vip_PE = c(
    "Slc1a3", "Riiad1", "Fads2", "Smim15", "Slc26a5", "Zc3hc1", 
    "Ppp1r1c", "Fsip1", "Prrg4", "Ptrh1", "Cd82", "Dffb", "Agtrap", 
    "Phf13", "Fam155b", "Slc35f6", "Dbx2", "Ccdc184", "Cgnl1", 
    "Plxnb1", "Gfral", "Esyt3", "Pml", "Rec114", "Ackr3", "Ing5", 
    "Sult1c2", "Tnfrsf12a", "Nat9", "Glp2r", "Aspm", "Slc25a37", 
    "Abhd4", "Pak1ip1", "Cebpd", "Insl3", "Lpl", "Prrc1", "Mmp2", 
    "Stox1", "Pop4", "Nsun2", "Rtl6", "Mbd3", "Tatdn2", "Pheta1", 
    "Nanp", "Znrd2", "Polr3g", "Fbxo32", "Dleu7", "Fos"
  ),
  Vip_PR = c(
    "Clcnkb", "RGD1310852", "Cd55", "Lpar6", "RGD1564786", 
    "Ror2", "Ca14", "Pgrmc2", "Ndst1", "RGD1565622", "Fos"
  )
)

# Perform GO analysis for each group and save results
go_results <- lapply(names(all_DEG_lists), function(group) {
  perform_go_analysis(all_DEG_lists[[group]], group)
})

# Create a workbook and add each group's results to a separate sheet
wb <- createWorkbook()
for (i in seq_along(go_results)) {
  sheet_name <- names(all_DEG_lists)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, go_results[[i]])
}

# Save the workbook
saveWorkbook(wb, file = "BLA_GO_analysis_results_all_combined.xlsx", overwrite = TRUE)




# Define gene lists for PC group
all_DEG_lists <- list(
  VGLUT1a_UE = c(
    "Arc", "Ecrg4", "Nr4a2", "Egr4", "Nr4a1", "Lonrf3", "Slc13a4", "Pde10a", 
    "Egr2", "Dusp6", "Nlrp12", "Nr4a3", "Mas1", "Egr1", "Pde7b", "Cox7a1", 
    "Bdnf", "Camk4", "Iyd", "Apoc2", "Myct1", "Ggn", "Ins2", "Cox6a2", 
    "Nlrp5", "Lsp1", "AC127887.1", "Galntl5", "Klre1", "Ly49si1", "Oxtr", 
    "Astl", "Lrrc26", "Cd52", "Cdkn2c", "Rpl35a.1", "Krtcap3", "RGD1560314", 
    "Syce3", "Epyc", "Aqp2", "Rdh8", "Plscr2", "AC105645.5", "RGD1305464", 
    "Fam124b", "RGD1307182", "AC132020.1", "Mroh9", "Stfa2l1", "Ropn1", 
    "Spaca7", "Ptpn20", "Bst2", "Nt5dc2", "RT1-CE11", "Gadd45b", "Tll1", 
    "Cntn1", "Baz1a", "Rnf17", "Pmfbp1", "Gpsm3", "Pdp1", "R3hdm1", "Pde4b", 
    "Cpne8", "Itgb2", "Lingo1", "Dusp2", "Sorcs1", "Frmpd4", "Rapgef5", 
    "Grik2", "Elmo1", "Ntrk2", "Bmp4", "Ppargc1a", "Vgf", "Grm1", "Gpr158", 
    "Dpyd", "Trerf1", "Cdh13", "Spaca5", "Lbhd2", "Clasp2", "Klhl2", "Fos"
  ),
  
  VGLUT1a_PE = c(
    "Nr4a3", "Nr4a1", "Dusp6", "Egr1", "Egr2", "Egr4", "Ntrk2", "Nr4a2", 
    "Pde10a", "Pde7b", "Fam110d", "R3hdm1", "Tiparp", "Phf21b", "Bdnf", 
    "Tll1", "Brinp1", "Arc", "Homer1", "AC120310.1", "Rapgef5", "Cenpm", 
    "Prok2", "Tpte2", "Gadd45b", "Vgf", "Gpr39", "Lonrf3", "Ralgapa2", "Mas1", 
    "Elmo1", "Ptgs2", "Hectd2", "Zdbf2", "Gpr158", "Ppp2r5a", "Arhgap42", 
    "Cdh13", "Nptx2", "Sorcs1", "Htra4", "Ms4a2", "Tpbgl", "Adora3", "Hspb3", 
    "Lpar3", "Tusc1", "Cib4", "Abhd12b", "Olr1078", "Epor", "Gpr62", "Adgrf1", 
    "Scimp", "RGD1559482.1", "Nr5a2", "Lgals3", "Olr1646", "Chodl", 
    "Slc25a6", "Hsd17b2", "RT1-Bb", "Oas1e", "Pde4b", "Pam", "Hivep1", 
    "Grm1", "Alkal2", "Npas4", "Osbpl8", "Trim9", "Numb", "Cnga4", "Plcl1", 
    "Slc24a3", "Spred3", "Kcnf1", "Mapk4", "Snai3", "Maml3", "Igf1r", 
    "Camk4", "Prickle1", "Tef", "Diras2", "Rell2", "Nfil3", "Prkca", 
    "Cckbr", "Fos"
  ),
  
  VGLUT1a_UR = c(
    "Npas4", "Nr4a1", "Egr4", "Dusp6", "Arc", "Vgf", "Egr1", "Zc3hav1", 
    "Zdbf2", "Ntrk2", "Etv5", "Nr4a3", "Clic3", "Col7a1", "Skida1", 
    "Pde10a", "Slc6a17", "Hmgcr", "Kcnj12", "Ncapg", "Tiparp", "R3hdm1", 
    "Il21r", "Acp7", "Cuzd1", "Fpr3", "Ffar2", "Rhod", "Olr210", "Slc44a3", 
    "Crhbp", "Cav2", "En2", "Gprc5a", "Galnt5", "Tktl1", "Radx", "Abcg8", 
    "RGD1305207", "Dnajc22", "Itga7", "Slurp1", "Cd3e", "Plet1", "RGD1562811", 
    "Wnt6", "Tbx2", "P2rx5", "Gast", "Steap3", "Lamb3", "Lrrc52", "Myoc", 
    "Tmprss11d", "Proser2", "Ppp1r3g", "Olr1568", "Olig2", "Lrp1b", "Fos"
  ),
  
  VGLUT1a_PR = c(
    "Egr2", "Arc", "Nr4a3", "Egr4", "Egr1", "Nr4a2", "Trib1", "Vgf", 
    "Gadd45g", "Gadd45b", "Hmga1", "Ntrk2", "Bdnf", "Homer1", "Dusp6", 
    "Nfil3", "Nr4a1", "Lingo1", "Alkal2", "Mas1", "R3hdm1", "Elmo1", 
    "Pde10a", "Cnksr1", "Idi1", "Irs2", "Abra", "Arl4d", "Sik1", "Spry2", 
    "Scg2", "Cdkn1a", "Tiparp", "Pim3", "Lonrf3", "Trim30c", "Cox6a2", 
    "Ccn1", "Klra17", "Tp53tg5", "Tfap2c", "RGD1564400", "Fam221b", 
    "Ly6c", "Ky", "Adgrf5", "RGD1559667", "Krt15", "Gast", "Hhipl2", 
    "Atf3", "Serpinb10", "Gnrhr", "Gch1", "Prl4a1", "RGD1561671", "Rab20", 
    "Msh5", "Sfta2", "Osbpl8", "Prx", "Actg2", "Cd200r1", "Nap1l5", 
    "Slc6a17", "Sorcs1", "Ly6k", "Amer3", "Arpp19.1", "Rapgef5", "Dnajb5", 
    "Numb", "Smad7", "Pde7b", "Atp1b1", "Trerf1", "Sik2", "Ryr2", "Cmtm5", 
    "Rheb", "Nrn1", "Grik2", "Hspa8", "Ralgapa2", "Lyz2", "Plk2", "Pank1", 
    "Csdc2", "Atp1a1", "Kcnv1", "Frmpd4", "Hmgcr", "Gpr158", "Pam", "Pim1", 
    "Frmd4a", "Cbarp", "Tuba1b", "Ppp6r2", "Lrrc8c", "Dnajc1", "Baz1a", 
    "Ier5", "Ptgs2", "Cplane1", "H3f3a", "Rps6ka2", "Fos"
  ),
  VGLUT1a_UR_no_common = c(
    "Npas4",  "Zc3hav1", 
    "Zdbf2",  "Etv5", "Clic3", "Col7a1", "Skida1", 
     "Slc6a17", "Hmgcr", "Kcnj12", "Ncapg", "Tiparp", 
    "Il21r", "Acp7", "Cuzd1", "Fpr3", "Ffar2", "Rhod", "Olr210", "Slc44a3", 
    "Crhbp", "Cav2", "En2", "Gprc5a", "Galnt5", "Tktl1", "Radx", "Abcg8", 
    "RGD1305207", "Dnajc22", "Itga7", "Slurp1", "Cd3e", "Plet1", "RGD1562811", 
    "Wnt6", "Tbx2", "P2rx5", "Gast", "Steap3", "Lamb3", "Lrrc52", "Myoc", 
    "Tmprss11d", "Proser2", "Ppp1r3g", "Olr1568", "Olig2", "Lrp1b"
  ),
  
  VGLUT1a_PR_no_common = c(
    "Egr2",  "Nr4a2", "Trib1",  
    "Gadd45g", "Gadd45b", "Hmga1", "Bdnf", "Homer1",  
    "Nfil3", "Lingo1", "Alkal2", "Mas1", "Elmo1", 
     "Cnksr1", "Idi1", "Irs2", "Abra", "Arl4d", "Sik1", "Spry2", 
    "Scg2", "Cdkn1a", "Tiparp", "Pim3", "Lonrf3", "Trim30c", "Cox6a2", 
    "Ccn1", "Klra17", "Tp53tg5", "Tfap2c", "RGD1564400", "Fam221b", 
    "Ly6c", "Ky", "Adgrf5", "RGD1559667", "Krt15", "Gast", "Hhipl2", 
    "Atf3", "Serpinb10", "Gnrhr", "Gch1", "Prl4a1", "RGD1561671", "Rab20", 
    "Msh5", "Sfta2", "Osbpl8", "Prx", "Actg2", "Cd200r1", "Nap1l5", 
    "Slc6a17", "Sorcs1", "Ly6k", "Amer3", "Arpp19.1", "Rapgef5", "Dnajb5", 
    "Numb", "Smad7", "Pde7b", "Atp1b1", "Trerf1", "Sik2", "Ryr2", "Cmtm5", 
    "Rheb", "Nrn1", "Grik2", "Hspa8", "Ralgapa2", "Lyz2", "Plk2", "Pank1", 
    "Csdc2", "Atp1a1", "Kcnv1", "Frmpd4", "Hmgcr", "Gpr158", "Pam", "Pim1", 
    "Frmd4a", "Cbarp", "Tuba1b", "Ppp6r2", "Lrrc8c", "Dnajc1", "Baz1a", 
    "Ier5", "Ptgs2", "Cplane1", "H3f3a", "Rps6ka2"
  ),
  
  VGLUT1b_PE = c(
    "Vom2r6", "Adgrg6", "Alpk3", "Ms4a4a", "Arhgap19", "Insc", "Smco2", 
    "Ermn", "Prr5l", "Lin28a", "Lpar1", "Cyp7a1", "Rs1", "RGD1308750", 
    "Papln", "AC132539.2", "Rhag", "Spata20", "Dppa1", "Stc2", "Foxc2", 
    "Irf2bpl", "Tiparp", "Egr2", "Inha", "Zbtb9", "Aldh16a1", "Med9", 
    "Papss2", "Spata13", "Nmbr", "Lrrc70", "Wwtr1", "Ccdc141", "Slc17a9", 
    "F2", "Lbh", "AC108595.1", "Stat6", "Ndrg1", "Arc", "Gfral", "Des", 
    "Trpm8", "Tekt4", "P3h4", "RGD1309106", "Scarf2", "Sh3bp4", "Rfesd", 
    "Masp2", "Tlcd3a", "Mas1", "Pdcl3", "Slc6a5", "Htra1", "Kash5", 
    "Ppp1r14a", "Fads2", "Cd53", "AC112531.1", "Hapln2", "Zc3hav1", 
    "Creb5", "Bcas1", "Usp51", "Tmprss6", "Oc90", "Ccnb2", "Zfp174", 
    "Grin2c", "Khnyn", "Gpr15", "Dusp4", "Spata4", "Slc29a3", "Ephb4", 
    "Adi1", "Lyrm2", "Nkain4", "Rad51ap1", "Tjap1", "Kcnj14", "Chek1", 
    "Fbxl21", "Plat", "Slitrk1", "Gtf3c4", "Fos"
  ),
  
  SL_UE = c("RGD1311946", "AC113910.1", "Orc1", "Zfp7", "Fos"),
  
  SL_UR = c("Fgd3", "Plekhb1", "Marveld1", "Jmjd4", "Pla2g1b", "Utp11", "Fos"),
  
  Vip_UE = c(
    "Rrad", "Ftsj3", "Cat", "Cdkn2aip", "Col4a3", "Pars2", "Adat3", 
    "Bmp1", "Rtl6", "Nfkbid", "Cdpf1", "Arrdc2", "Pgap2", "Imp4", "Nr4a1", 
    "Aldh3b2", "Bag3", "Slc22a12", "AC107531.3", "Slc1a5", "Plg", "Fermt3", 
    "Rab3il1", "Clca5", "Mcm2", "Tprn", "Emilin3", "Ap5s1", "Ntng2", 
    "Olr633", "Cerkl", "Hmgb4", "Crispld1", "Cnksr1", "Aunip", "Ptgr1", 
    "Abcb5", "Gadd45b", "Aox4", "Tfeb", "Ppl", "Mgrn1", "Ccdc40", "Gemin4", 
    "Olr1389", "Phospho1", "Ccl5", "Areg", "Fam3d", "Piwil2", "Fgfr4", 
    "Ctla2a", "Akr1c1", "Chodl", "Slc25a15", "Afap1l1", "Pou4f3", 
    "Arl14epl", "Chst9", "Piezo1", "Irf8", "Olr1734", "Psmb9", "Tex26", 
    "Opa3", "Tmem47", "Septin1", "Adgrv1", "Nudt10", "Egr4", "Vstm4", 
    "Tpd52l1", "Fos"
  ),
  
  Vip_PE = c(
    "Vangl1", "Rps18l1", "Tedc2", "Bcl6", "Erf", "Clec1a", "Cyp1b1", 
    "Spc24", "Zp3r", "Ccn6", "Nr4a3", "Rgs19", "Gsta1", "Mon1b", "Avpi1", 
    "Tbl3", "Tlcd5", "Vom2r35", "Fpr2", "Cdhr5", "Asah2", "Kctd15", 
    "Olr324", "Wnt8b", "Gnat2", "Tnfsf10", "Tmem144", "Kcna5", "Arhgef16", 
    "Tas1r3", "Actl7a", "Hcrtr1", "Ca8", "Agmo", "Tjp3", "Bcl2a1", 
    "Snx33", "Pls1", "Cyp39a1", "AC132020.1", "Ccdc40", "Pdia2", "Elf3", 
    "Cfhr2", "Sele", "Myh6", "Ccdc74a", "Neil3", "Pde6a", "Psmb9", "Ggt1", 
    "Gtf3c6", "E2f2", "Pon2", "Ninl", "Glipr2", "Mob3c", "Hsd17b13", 
    "Mafk", "Ogfod1", "Mfsd5", "Sdhaf3", "Nacc2", "Pnpla1", "Fos"
  ),
  
  Vip_UR = c(
    "Egr4", "Lck", "Slc25a24", "Tekt2", "Mtmr10", "Alg10", "Slc31a2", 
    "Kcng2", "Urm1.1", "Wnt4", "Tnfsf12", "Prdm11", "Fos"
  ),
  
  Vip_PR = c("Dqx1", "Pars2", "Fos"),
  
  Sst_UE = c("Sco1", "Fos"),
  
  Sst_PE = c(
    "Ybey", "Zpbp2", "Fpr2", "RGD1310166", "Cd68", "AC130232.2", "Spatc1l", 
    "Fxyd6", "Tcf7l2", "Prok2", "Aurka", "Ell3", "Nr4a1", "Ccdc159", 
    "Cndp1", "Cnn2", "Fos"
  ),
  
  Sst_UR = c(
    "Rgs20", "Ogg1", "Gipr", "Fhl3", "Dpy30", "Wnt9b", "Colgalt1", 
    "Ghr", "Bean1", "Plce1", "Trim36", "Adra2a", "Hyal2", "Rlbp1", 
    "Tmem159", "Plaat5", "Zfp764l1", "Mvp", "Saxo2", "Ifitm6", "Fabp12", 
    "Crh", "Gpr160", "Si", "Mgst2", "Ret", "Alox5", "Cav1", "Tspan11", 
    "Abcc9", "Mgst1", "Bcas1", "Lrrc19", "Epas1", "Zdhhc22", "Meox2", 
    "Wnt7b", "Gpr182", "Phf21b", "Has2", "Notch3", "Kank2", "Drd2", 
    "Dpy19l2", "RGD1560917", "Abhd14b", "Spink8", "Otos", "Znrf4", 
    "Vom2r76", "Kif6", "Crhr1", "Tcam1", "Myo15a", "Gpr37l1", "Cdh19", 
    "Kcnj10", "C1qtnf7", "Spata13", "Cysltr2", "Foxc1", "Parp14", 
    "Thpol1", "Olig2", "Pou1f1", "Onecut2", "Dipk1c", "RT1-N3", "Ybey", 
    "Ocm", "Vom2r77", "Pcdhb22", "Dffb", "Fos"
  ),
  
  Sst_PR = c(
    "Slc66a3", "Ccne1", "Zcchc3", "Trim36", "Ptprh", "Ceacam1", "Fam181b", 
    "Cetn4", "Anxa4", "Olr854", "Plcd1", "Abca8a", "Hist1h2bo", "Wnt5a", 
    "Rrad", "Fos"
  ),
  
  IN1_PE = c(
    "Tmem190", "Hif3a", "Ccdc88b", "Bglap", "F2rl1", "Dkk2", "Them5", 
    "Tmem229a", "Gimap6", "Galnt5", "Nsun4", "Krt8", "Notch3", "Hs3st3b1", 
    "Itgb3", "Nr5a2", "Prl3d4", "Slc39a12", "Tex43", "Tcf19", "Cdc45", 
    "Tmem202", "Lrrc39", "Tpm2", "Atoh8", "Mdfic", "Trpv5", "Cdh26", 
    "Cdca8", "Lhcgr", "Hephl1", "Icam4", "Sox15", "Dnah14", "RGD1561440", 
    "P2rx6", "Emp3", "Ceacam9", "Map7d3", "Hmga2", "Ttc12", "Ggt6", 
    "Tmprss7", "F12", "S100a3", "Atic", "Yap1", "Mfap4", "Crhbp", "Dynlt2", 
    "Apom", "Ttll9", "Aldh3b2", "Gcnt4", "Lrrn4", "Depdc7", "Cyp2j16", 
    "Plcd4", "RGD1309106", "Pdha1l1", "Cngb1", "Acacb", "Rhcg", "Nrarp", 
    "Asap3", "Spag8", "Ccdc154", "Spta1", "Cd38", "Sh3tc1", "Hormad2", 
    "Thbs3", "Trpv6", "Pigo", "Trpc5os", "Bard1", "Prl5a1", "Agrp", 
    "RGD1305464", "Dkk4", "Tap2", "Gpt2", "Gjb6", "Pard6b", "Rfxap", 
    "Pecr", "Notch4", "Ctxnd1", "Etv5", "Hemk1", "Etv4", "Pcdhb20", 
    "Nxnl1", "Gsta6", "Serpinh1", "Mettl1", "Irf2bpl", "Pgghg", "Acss3", 
    "Wif1", "Cand2", "Shh", "Arhgef19", "Bmp4", "Pcolce", "Relt", "Mt1", 
    "Axl", "Caskin2", "RGD1564899", "Pop7", "Vwc2l", "Plpp3", "Iqsec3", 
    "Abcc8", "Helz", "Sema3e", "Tmem64", "Lats2", "Hapln4", "Kcnh8", 
    "Cptp", "Tpmt", "Med17", "Recql", "Pnpo", "Notch1", "Ifrd1", "Pygm", 
    "Psmd12", "Chaf1a", "Hunk", "Mtrf1l", "Slc16a6", "Pros1", "Plpbp", 
    "Fam162a", "Il6st", "Trpc3", "Map1s", "Slc6a9", "Wscd1", "Abcb9", 
    "Cd300lf", "Eif4ebp1", "Ifitm3", "Nr1h3", "Adamts15", "Vstm4", 
    "Nudt1", "Mmp23", "C1qtnf5", "Gfral", "Kdr", "Parp9", "Aqp1", "Ucn", 
    "Ttll8", "Stk17b", "F5", "Slc4a2", "RGD1560436", "Tigd3", "Neu2", 
    "Adap2", "Cdh19", "Lypd6", "Tbc1d8b", "Hba-a3", "Acadm", "Lrp11", 
    "Ccdc157", "Pias4", "Glis2", "Timm44", "Bcl6", "Nr2c1", "Syt3", 
    "Tmem184b", "Kcnv1", "Ube2e1", "Nxpe3", "Msantd4", "Chmp3", "Lamtor3", 
    "Castor2", "Dhtkd1", "Mpdu1", "Fos"
  ),
  
  IN1_UR = c(
    "Dusp5", "Zfp36", "Pyroxd2", "Esrra", "Npas4", "RGD1565057", "Wnt2b", 
    "Tacr3", "Kcnmb2", "AC111647.1", "Prok2", "Zfp786", "Traf2", "Insm1", 
    "Rassf2", "AC127963.1", "Arhgap36", "Arc", "Mns1", "Col6a6", "Trarg1", 
    "Pmp22", "Pdc", "Pogk", "Amtn", "Neil2", "Dok2", "RGD1309651", 
    "Col13a1", "Nudt1", "Ptgfrn", "Tmem251", "Gigyf1", "Ephx1", "Tec", 
    "Mthfd1l", "Msantd3", "Psmc3ip", "Nr4a1", "Ankdd1a", "Tek", "Adora2a", 
    "Ube2ql1", "Smad9", "Abcb4", "Sspn", "Kbtbd8", "Chadl", "Sntb1", 
    "Drd2", "Lca5", "Slc2a12", "Mturn", "Ly75", "Tmem51", "Oma1", 
    "Mfap4", "AC131411.2", "Haus4", "Gng4", "Htra4", "Erich1", "Cfap53", 
    "Pdp2", "Med20", "Mrpl19", "Spdya", "Iah1", "Med25", "Cldn12", 
    "St8sia3", "Syt12", "Pak2", "Atox1", "N4bp1", "Spsb3", "Dhrs7", 
    "Penk", "Scg2", "Hpca", "Gnl1", "Dnajc14", "Pitrm1", "Guf1", "Rgs2", 
    "Rhbdd2", "Hdx", "Amt", "Polk", "Cdca4", "Osbp", "Pramef8", "Lmod1", 
    "Adamts3", "Dhodh", "Lrrc43", "Ankrd46", "Mfsd10", "Dcbld2", 
    "Agtrap", "Nts", "Marchf10", "Ddx49", "Hps6", "Nkpd1", "Adra1d", 
    "Ak8", "Mid1", "Slco2a1", "Krt26", "Itgae", "Cd55", "Adamts14", 
    "Uros", "Peli3", "Dolpp1", "Fblim1", "Csrnp1", "Ebi3", "Esrrg", 
    "Diaph3", "Pwwp2b", "Pdzk1", "Cacna2d4", "Arx", "Hyal2", "Btbd17", 
    "Heatr6", "Zfp951", "Cln8", "Fbn2", "Snu13.1", "Caskin1", "Hfm1", 
    "Speg", "Fance", "Lpcat1", "Dhx35", "Ctsb", "Ttll4", "Kcnab1", 
    "Arf5", "Gnb5", "Tmem30a", "Rasl11b", "Efna3", "Actn2", "Acot9", 
    "Nup205", "Kdm1a", "Grpel1", "Gpr19", "Mrpl32", "Mitd1", "Tbc1d30", 
    "Tmem167b", "Plekhh1", "Il16", "Fos"
  ),
  
  IN2_UE = c(
    "Col16a1", "Selenov", "Ldhc", "Vwce", "Bcam", "Cep55", "Zfp663", 
    "Sun5", "Actl7b", "Grpr", "Rasl12", "Gpr62", "Ets1", "Chrna3", 
    "Cdrt4", "C1qtnf1", "Gdf9", "Tekt3", "Slc16a3", "Fam72a", "Ddr2", 
    "Niban1", "Phf11", "Slc17a2", "Armc3", "Map3k8", "Stfa2l1", "Slc51a", 
    "Htr1f", "Calhm4", "Ubash3a", "Mettl18", "Vsir", "Gpr37", "Rp2", 
    "Ube2g2", "Wdr88", "Aspdh", "Abca1", "Socs1", "Lrrc46", "Dkk2", 
    "Gsta4", "Zswim7", "Rhoj", "Pappa2", "Trim15", "Fos"
  ),
  
  IN2_PE = c(
    "Egr1", "Masp1", "Klf10", "Ptprcap", "RGD1560303", "Cep55", "Pdilt", 
    "Magel2", "Dnase2b", "Dusp14l1", "Ret", "Col1a2", "Serinc4", "Emilin3", 
    "Awat1", "Fshr", "Il2rb", "Angptl4", "Cyth4", "Cd3g", "Treml1", 
    "Abca6", "Myocd", "Adora2b", "Stc2", "Ube2d4", "Rtp4", "Slc25a15", 
    "Proc", "Ces4a", "Crybb1", "Tctn3", "Myo1e", "Nr4a1", "Itprip", 
    "Tbxa2r", "N4bp3", "Gmppb", "Lbh", "Inpp5j", "Fos"
  ),
  
  IN2_PR = c(
    "Igfbpl1", "Etv4", "Cspg4b", "Sema5b", "Npas4", "Ubtd1", "Il31ra", 
    "Mttp", "Pon2", "Lmcd1", "Sbspon", "Acnat1.1", "Cyld-ps1", "Zfp449", 
    "Prox2", "Gpr182", "AC096792.1", "Gtf2a2", "Parp3", "Cx3cr1", 
    "Vom2r77", "Wnt6", "Ccr10", "Avpr1b", "Fcgr3a", "Pkd1l1", "Cby2", 
    "Nkapl", "Gpr15", "Fbn2", "Slc22a16", "Fos"
  )
)

# Perform GO analysis for each group and save results
go_results <- lapply(names(all_DEG_lists), function(group) {
  perform_go_analysis(all_DEG_lists[[group]], group)
})

# Create a workbook and add each group's results to a separate sheet
wb <- createWorkbook()
for (i in seq_along(go_results)) {
  sheet_name <- names(all_DEG_lists)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, go_results[[i]])
}

# Save the workbook
saveWorkbook(wb, file = "PC_GO_analysis_results_all_combined.xlsx", overwrite = TRUE)



# Read the Excel sheets into R
file_path <- "pathyway/BLA UR PR GO.xlsx"
VGLUT1_UR <- read_excel(file_path, sheet = "VGLUT1_UR")
VGLUT1_PR <- read_excel(file_path, sheet = "VGLUT1_PR")

# Add a 'Condition' column to each dataframe
VGLUT1_UR$Condition <- "VGLUT1_UR"
VGLUT1_PR$Condition <- "VGLUT1_PR"



# Combine the datasets (assuming you've already filtered for red-highlighted GO terms)
combined_data <- rbind(VGLUT1_UR,VGLUT1_PR)

# Convert the 'Condition' column to a factor with specific levels
combined_data$Condition <- factor(combined_data$Condition, levels = c("VGLUT1_UR", "VGLUT1_PR"))


# Create the bubble plot with swapped legend positions and fewer size scale breaks
VGLUT_bubble <- ggplot(combined_data, aes(x = Condition, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(6, 12), name = "Count", breaks = c(2, 4, 6)) +  # Fewer breaks for size legend
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(title = "VGLUT1 significant GO Terms",
       x = "Condition",
       y = "GO Terms",
       color = "Adjusted P-value",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")) +
  guides(color = guide_colorbar(order = 1), # Adjusted p-value legend first
         size = guide_legend(order = 2))    # Count legend second

# Save the plot with adjusted legends and fewer size breaks
ggsave("BLA_VGLUT1_Bubble v2.tiff", plot = VGLUT_bubble, width = 10, height = 10, dpi = 300, units = "in", bg = "white")




# Read the Excel sheets into R
file_path <- "pathyway/BLA bubble2.xlsx"

VGLUT1_PR <- read_excel(file_path, sheet = "VGLUT1_PR")
Sst_PR <- read_excel(file_path, sheet = "Sst_PR")
Npy_PR <- read_excel(file_path, sheet = "Npy_PR")



# Add a 'Condition' column to each dataframe

VGLUT1_PR$Condition <- "VGLUT1_PR"
Sst_PR$Condition <- "Sst_PR"
Npy_PR$Condition <- "Npy_PR"



# Combine the datasets (assuming you've already filtered for red-highlighted GO terms)
combined_data <- rbind(Sst_PR,Npy_PR,VGLUT1_PR)

# Convert the 'Condition' column to a factor with specific levels
combined_data$Condition <- factor(combined_data$Condition, levels = c("Npy_PR","Sst_PR","VGLUT1_PR"))


# Create a blank row with NA values for all columns
blank_row <- combined_data[1, ]
blank_row[,] <- NA  # Set all values to NA

# Fill the 'Condition' and 'Description' columns appropriately
blank_row$Condition <- "VGLUT1_PR" # or any Condition level
blank_row$Description <- ""  # or leave it as NA

# Append the blank row to the end of the combined data
combined_data <- rbind(combined_data, blank_row)


# Create the bubble plot for inhibitory neurons
VGLUT_bubble <- ggplot(combined_data, aes(x = Condition, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(6, 12), name = "Count", breaks = c(2, 4, 6)) +  # Fewer breaks for size legend
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(title = "VGLUT1 significant GO Terms",
       x = "Condition",
       y = "GO Terms",
       color = "Adjusted P-value",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")) +
  guides(color = guide_colorbar(order = 1), # Adjusted p-value legend first
         size = guide_legend(order = 2))    # Count legend second

# Save the plot with adjusted legends and fewer size breaks
ggsave("BLA VGLUT1 Sst Npy bubble.tiff", plot = VGLUT_bubble, width = 10, height = 10, dpi = 300, units = "in", bg = "white")





# Read the Excel sheets into R
file_path <- "pathwayPC_bubble.xlsx"

UE <- read_excel(file_path, sheet = "VGLUT1a_UE")
PE <- read_excel(file_path, sheet = "VGLUT1a_PE")
UR <- read_excel(file_path, sheet = "VGLUT1a_UR")
PR <- read_excel(file_path, sheet = "VGLUT1a_PR")

# Add a 'Condition' column to each dataframe
UE$Condition <- "VGLUT1a_UE"
PE$Condition <- "VGLUT1a_PE"
UR$Condition <- "VGLUT1a_UR"
PR$Condition <- "VGLUT1a_PR"



# Create a blank row with NA values for all columns
blank_row <- combined_data[1, ]
blank_row[,] <- NA  # Set all values to NA

# Fill the 'Condition' and 'Description' columns appropriately
blank_row$Condition <- "VGLUT1a_UE" # or any Condition level
blank_row$Description <- ""  # or leave it as NA


# Combine the datasets (assuming you've already filtered for red-highlighted GO terms)
combined_data <- rbind(UE,PE,UR,PR)

# Convert the 'Condition' column to a factor with specific levels
combined_data$Condition <- factor(combined_data$Condition, levels = c("VGLUT1a_UE", "VGLUT1a_PE", "VGLUT1a_UR", "VGLUT1a_PR"))


# Create the bubble plot
VGLUT_bubble <- ggplot(combined_data, aes(x = Condition, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10), name = "Count") +
  scale_color_gradient(low = "red", high = "blue",name = "Adjusted p-value") +
  labs(title = "PC pyramidal significant GO Terms",
       x = "Condition",
       y = "GO Terms",
       color = "Adjusted P-value",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  guides(color = guide_colorbar(order = 1), # Adjusted p-value legend first
       size = guide_legend(order = 2))    # Count legend second

ggsave("PC_pyramidal_Bubble v3.tiff", plot  = VGLUT_bubble, width = 10, height = 10, dpi = 300, units = "in",bg = "white")





# BLA, PC comparison, Read the Excel sheets into R
file_path <- "pathyway/BLA vs PC excitatory bubble.xlsx"

BLA_PR <- read_excel(file_path, sheet = "VGLUT1_PR")

PC_PR <- read_excel(file_path, sheet = "VGLUT1a_PR")

# Add a 'Condition' column to each dataframe

BLA_PR$Condition <- "VGLUT1_PR"

PC_PR$Condition <- "VGLUT1a_PR"




# Combine the datasets (assuming you've already filtered for red-highlighted GO terms)
combined_data <- rbind(BLA_PR,PC_PR)

# Convert the 'Condition' column to a factor with specific levels
combined_data$Condition <- factor(combined_data$Condition, levels = c( "VGLUT1_PR",  "VGLUT1a_PR"))


# Create the bubble plot
VGLUT_bubble <- ggplot(combined_data, aes(x = Condition, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10), name = "Count") +
  scale_color_gradient(low = "red", high = "blue",name = "Adjusted p-value") +
  labs(title = "BLA versus PC  GO Terms",
       x = "Condition",
       y = "GO Terms",
       color = "Adjusted P-value",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  guides(color = guide_colorbar(order = 2), # Adjusted p-value legend first
       size = guide_legend(order = 1))    # Count legend second

ggsave("BLA_PC_comparison_Bubble v2.tiff", plot  = VGLUT_bubble, width = 10, height = 10, dpi = 300, units = "in",bg = "white")




