



#### COST PARTITION CALCULATION ####
# 
# 
# 
# cost_df <- df_param_kmax_alpha %>% 
#   left_join(df %>% 
#               dplyr::select(Species, source, Drydown.days,ca,T,D,Iabs_used,gC) %>% 
#               group_by(Species, source) %>% 
#               summarise_all(mean, na.rm=TRUE)) %>% 
#   rowwise() %>%
#   do(mapped = get_cost(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
#   unnest(cols = c(mapped))
#   
# 
# 
# cost_df %>% rename(psi_s = var) %>%  left_join(df_param_kmax_alpha,by=c("Species","source","scheme")) %>% 
#   group_by(Species,scheme,source) %>%  
#   mutate(h_cost = case_when(scheme=="SOX"~1-h_cost/a,
#                             TRUE~h_cost),
#          total_cost = jmax_cost + h_cost,
#          jmax_cost = jmax_cost/max(jmax_cost,na.rm=TRUE),
#          h_cost = h_cost/max(h_cost,na.rm=TRUE)) %>% 
#   ggplot(aes(psi_s,(jmax_cost-h_cost),group=interaction(Species,source,scheme)))+ 
#   geom_line(color="grey20") +
#   geom_hline(yintercept = 0, linetype=2,color="grey50")+
#   facet_grid(Species~scheme,scales ="free")+
#   theme_bw()
# 
# cost_df %>%
#   group_by(Species,scheme,source) %>%
#   mutate(h_cost = case_when(scheme=="SOX"~1-h_cost/a,
#                             scheme=="PROFITMAX2"~h_cost/a,
#                             TRUE~h_cost),
#          total_cost = jmax_cost + h_cost,
#          jmax_cost = jmax_cost/max(jmax_cost,na.rm=TRUE),
#          h_cost =h_cost/max(h_cost,na.rm=TRUE)) %>%
#   ggplot(aes(group=interaction(Species,source,scheme)))+
#   geom_line(aes(var,h_cost)) +
#   geom_line(aes(var,jmax_cost), color="red") +
#   # ylim(0,100)+
#   facet_grid(Species~scheme,scales ="free")+
#   theme_bw()
# # 
# # cost_df %>% 
# #   group_by(Species,scheme,source) %>% 
# #   mutate(h_cost = case_when(scheme=="SOX"~h_cost/a,
# #                             TRUE~h_cost),
# #          total_cost = jmax_cost + h_cost,
# #          jmax_cost = jmax_cost/max(total_cost,na.rm=TRUE),
# #          h_cost =h_cost/max(total_cost,na.rm=TRUE)) %>% 
# #   ggplot(aes(group=interaction(Species,source,scheme)))+ 
# #   geom_line(aes(var,h_cost)) +
# #   geom_line(aes(var,jmax_cost), color="red") +
# #   # ylim(0,100)+
# #   facet_grid(Species~scheme,scales ="free")
# 
# 
# cost_df %>% rename(psi_s = var) %>%  left_join(df_param_kmax_alpha,by=c("Species","source","scheme")) %>% 
#   group_by(Species,scheme,source) %>% 
#   mutate(h_cost = case_when(scheme=="SOX"~h_cost/a,
#                             TRUE~h_cost),
#          total_cost = jmax_cost + h_cost,
#          jmax_cost = jmax_cost/total_cost,
#          h_cost =h_cost/total_cost) %>% 
#   ggplot(aes(group=interaction(Species,source,scheme)))+ 
#   geom_line(aes(psi_s,h_cost)) +
#   geom_line(aes(psi_s,jmax_cost), color="red") +
#   geom_vline(aes(xintercept = P50))+
#   # ylim(0,100)+
#   facet_grid(Species~scheme,scales ="free")+
#   theme_bw()
#   



#### SIMULATION + PARAMETERS ANALYSIS ####

# df_a_diff <- df %>% 
#   group_by(scheme,acclimation,Species,source) %>% 
#   filter(!is.na(A)) %>% 
#   dplyr::select(-c(dpsi,genus,species,subsp,P88..MPa.,P50..MPa.,P12..MPa.,
#              SLA..cm2.g.1.,Height.max..m.,Ks..kg.m.1.MPa.1.s.1.,Huber.value,
#              KL..kg.m.1.MPa.1.s.1.)) %>% 
#   mutate(diff_a = A - a_pred) %>% 
#   summarise(n_dist = n(),
#             r = cor(A, a_pred, use = "pairwise.complete.obs"),
#             bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
#             rmse = Metrics::rmse(A,a_pred),
#             beta = lm(A~a_pred)$coefficients[2]) %>%
#   left_join(df_param%>%
#               dplyr::select(- c(dpsi, inst,genus,species,subsp,P88..MPa.,P50..MPa.,P12..MPa.,
#                          SLA..cm2.g.1.,Height.max..m.,Ks..kg.m.1.MPa.1.s.1.,Huber.value,
#                          KL..kg.m.1.MPa.1.s.1.))) %>%
#   ungroup() %>% 
#   pivot_wider(names_from = acclimation, values_from = c(6:14)) %>% 
#   mutate(r_diff = r_Acclimated - `r_No acclimated`,
#          bias_diff = bias_Acclimated - `bias_No acclimated`,
#          rmse_diff = rmse_Acclimated - `rmse_No acclimated`,
#          beta_diff = beta_Acclimated - `beta_No acclimated`,
#          K.scale_diff = K.scale_Acclimated - `K.scale_No acclimated`,
#          psi50_diff = P50_Acclimated - `P50_No acclimated`,
#          b_diff = b_Acclimated - `b_No acclimated`,
#          gamma_diff = gamma_Acclimated - `gamma_No acclimated`
#   ) 
# 
# # A - r Pearson's
# p1 <- df_a_diff %>% 
#   # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
#   mutate(r_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~r_diff,
#                             TRUE~NA_real_)) %>% 
#   ggplot(aes(K.scale_diff, r_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*"K"))
# 
# p2 <- df_a_diff %>% 
#   # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
#   ggplot(aes(gamma_diff, r_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*gamma))
# 
# 
# annotate_figure(
#   ggarrange(p1,p2, 
#           align='hv', labels=c('a', 'b'),
#           common.legend = T,ncol=2, nrow = 1),
#   left = textGrob("(Acclimated - No acclimated) A rate r Pearson's correlation", 
#                   rot = 90, vjust = 1, gp = gpar(cex = 1)))
# 
# # A - BETA
# p1 <- df_a_diff %>% 
#   # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
#   mutate(beta_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~beta_diff,
#                                TRUE~NA_real_)) %>% 
#   ggplot(aes(K.scale_diff, beta_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*"K"))
# 
# p2 <- df_a_diff %>% 
#   # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
#   ggplot(aes(gamma_diff, beta_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*gamma))
# 
# annotate_figure(
#   ggarrange(p1,p2,
#             align='hv', labels=c('a', 'b'),
#             common.legend = T,ncol=2, nrow = 1),
#   left = textGrob("(Acclimated - No acclimated) A rate \u03B2", 
#                   rot = 90, vjust = 1, gp = gpar(cex = 1)))
# 
# 
# # A - RMSE
# p1 <- df_a_diff %>% 
#   # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
#   mutate(rmse_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~rmse_diff,
#                             TRUE~NA_real_)) %>% 
#   ggplot(aes(K.scale_diff, rmse_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*"K"))
# 
# p2 <- df_a_diff %>% 
#   # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
#   ggplot(aes(gamma_diff, rmse_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*gamma))
# 
# annotate_figure(
#   ggarrange(p1,p2,
#             align='hv', labels=c('a', 'b'),
#             common.legend = T,ncol=2, nrow = 1),
#   left = textGrob("(Acclimated - No acclimated) A rate RMSE", 
#                   rot = 90, vjust = 1, gp = gpar(cex = 1)))
# 
# 
# # A - BIAS
# p1 <- df_a_diff %>% 
#   # filter(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")) %>%
#   mutate(bias_diff = case_when(scheme %in% c("PROFITMAX","PROFITMAX2","SOX")~bias_diff,
#                                TRUE~NA_real_)) %>% 
#   ggplot(aes(K.scale_diff, bias_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*"K"))
# 
# p2 <- df_a_diff %>% 
#   # filter(scheme %in% c("WUE","CGAIN","CMAX",'PHYDRO')) %>% 
#   ggplot(aes(gamma_diff, bias_diff, color = scheme, group = scheme))+
#   geom_hline(yintercept = 0, color = "grey20", linetype = 3)+
#   geom_vline(xintercept = 0, color = "grey20", linetype = 3)+
#   # geom_encircle()+
#   geom_convexhull(alpha=0, fill = "transparent", size = 1.3, show.legend = FALSE)+
#   geom_point()+
#   # geom_density_2d()+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   # scale_color_manual(values = brewer_pal(palette = "RdYlGn")(4))+
#   # ylab(expression("(Acclimated - No acclimated) A rate r Pearson's correlation"))+
#   ylab("")+
#   xlab(expression(Delta*gamma))
# 
# annotate_figure(
#   ggarrange(p1,p2,
#             align='hv', labels=c('a', 'b'),
#             common.legend = T,ncol=2, nrow = 1),
#   left = textGrob("(Acclimated - No acclimated) A rate BIAS", 
#                   rot = 90, vjust = 1, gp = gpar(cex = 1)))
# 
# 
# 



################################################################################

# 
# Qi_params_no <- df_a_fix_param%>% 
#   filter(Species=="Quercus ilex",
#          acclimation == "No acclimated",
#          source == "Epron and Dreyer (1990)",
#          dpsi == FALSE)
# 
# Ep_params_no <- df_a_fix_param%>% 
#   filter(Species=="Eucalyptus pilularis",
#          acclimation == "No acclimated",
#          dpsi == FALSE)
# 
# Epo_params_no <- df_a_fix_param%>% 
#   filter(Species=="Eucalyptus populnea",
#          acclimation == "No acclimated",
#          dpsi == FALSE)
# 
# 
# df_sim_Qi_no <- Qi_params_no %>%
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)
# 
# df_sim_Ep_no <- Ep_params_no %>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)
# 
# df_sim_Epo_no <- Epo_params_no %>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_no_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)
# 
# df_summary_no <- df_sim_Qi_no  %>% 
#   rbind(df_sim_Ep_no) %>% 
#   rbind(df_sim_Epo_no) %>% 
#   # filter(a>0.2) %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)


# 
# Q_params <- df_a_fix_param%>% 
#   filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
#          acclimation == "Acclimated",
#          source == "Epron and Dreyer (1990)",
#          dpsi == FALSE)
# 
# Q_params_actual <- df_a_fix_param%>% 
#   left_join(df %>% 
#               dplyr::select(Species,ca,T,D,Iabs_used,gC) %>% 
#               group_by(Species) %>% 
#               summarise_all(mean, na.rm=TRUE)) %>% 
#   filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
#          acclimation == "Acclimated",
#          source == "Epron and Dreyer (1990)",
#          dpsi == FALSE)
# 
# Q_params_actual_kmax_alpha <- df_param_kmax_alpha %>% 
#   left_join(df %>% 
#               dplyr::select(Species,ca,T,D,Iabs_used,gC) %>% 
#               group_by(Species) %>% 
#               summarise_all(mean, na.rm=TRUE)) %>% 
#   filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
#          acclimation == "Acclimated",
#          source == "Epron and Dreyer (1990)",
#          dpsi == FALSE)
# 
# Ep_params <- df_a_fix_param%>% 
#   filter(Species=="Eucalyptus pilularis",
#          acclimation == "Acclimated",
#          dpsi == FALSE)
# 
# Epo_params <- df_a_fix_param%>% 
#   filter(Species=="Eucalyptus populnea",
#          acclimation == "Acclimated",
#          dpsi == FALSE)


# df_sim_Qi <- Qi_params %>%
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)
# 
# df_sim_Ep <- Ep_params %>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)
# 
# df_sim_Epo <- Epo_params%>% 
#   split(seq(nrow(.))) %>% 
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 4000)

# df_test <- df %>% 
#   filter(Species%in%c("Eucalyptus populnea",
#                       "Eucalyptus pilularis",
#                       "Quercus ilex"))
# 
# df_summary <- df_sim_Qi  %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   # filter(a>0.2) %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)

# 
# df_sim_Q_low <- Q_params %>%
#   split(seq(nrow(.))) %>%
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=1000,co2_now = 400)
# df_summary_low <- df_sim_Q_low  %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# df_sim_Q_high <- Q_params %>%
#   split(seq(nrow(.))) %>%
#   purrr::map_df(get_acclimated_response,tc_now=25,ppfd_now=1200,vpd_now=3000,co2_now = 400)
# df_summary_high <- df_sim_Q_high  %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# df_sim_Q_actual <- Q_params_actual %>%
#   rowwise() %>%
#   do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
#   unnest(cols = c(mapped))
# df_summary_actual<- df_sim_Q_actual  %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# df_sim_Q_actual_kmax_alpha <- Q_params_actual_kmax_alpha %>%
#   rowwise() %>%
#   do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
#   unnest(cols = c(mapped))
# df_summary_actual_kmax_alpha <- df_sim_Q_actual_kmax_alpha  %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# df_test <- df %>% 
#   filter(Species%in%c("Quercus ilex","Quercus petraea","Quercus pubescens"),
#          acclimation == "Acclimated",
#          source == "Epron and Dreyer (1990)")
# df_summary_sp <- df_test %>% 
#   group_by(Species) %>% 
#   dplyr::select(P50) %>% 
#   summarise_all(mean,na.rm=TRUE)


# 
# p1 <- df_sim_Q_low  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,jmax,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #            data=df_summary_low,
#   #            linetype = 1,size=1,
#   #            show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #            linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                     values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,170)+
#   NULL
# 
# p2 <- df_sim_Q_high  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,jmax,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary_high,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,170)+
#   NULL
# 
# p3 <- df_sim_Q_actual  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,jmax,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary_high,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,60)+
#   NULL
# 
# p4 <- df_sim_Q_actual_kmax_alpha %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,jmax,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 190, yend= 170,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary_high,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=190,yend=170),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,70)+
#   NULL
# 
# ggarrange(p1,p2, 
#           align='hv', labels=c('a', 'b'),
#           common.legend = T,ncol=2, nrow = 1)
# 
# 
# 
# p1 <- df_sim_Q_low  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,gs,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
#                arrow = arrow(length = unit(0.3, "cm")),
#                lineend = c('round'),linejoin = c('mitre'),
#                data=df_summary_low,
#                linetype = 1,size=1,
#                show.legend = FALSE)+
#   geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
#                arrow = arrow(length = unit(0.3, "cm")),
#                lineend = c('round'),linejoin = c('mitre'),
#                linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,0.6)+
#   NULL
# 
# p2 <- df_sim_Q_high  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,gs,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
#                arrow = arrow(length = unit(0.3, "cm")),
#                lineend = c('round'),linejoin = c('mitre'),
#                data=df_summary_high,
#                linetype = 1,size=1,
#                show.legend = FALSE)+
#   geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
#                arrow = arrow(length = unit(0.3, "cm")),
#                lineend = c('round'),linejoin = c('mitre'),
#                linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   ylim(0,0.6)+
#   NULL
# 
# p3 <- df_sim_Q_actual  %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,gs,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,gs,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1,alpha=0.5)+
#   geom_point(aes(LWP,gC),data=df_test)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 0.24, yend= 0.22,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary_high,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.24,yend=0.22),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   xlim(-6,0)+
#   # ylim(0,0.31)+
#   NULL
# 
# ggarrange(p1,p2, 
#           align='hv', labels=c('a', 'b'),
#           common.legend = T,ncol=2, nrow = 1)
# 
# p_dpsi_accl <- df_sim_Qi %>%
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,dpsi,color = scheme),size = 1)+
#   # geom_line(aes(p_leaf,dpsi,color = scheme,group = interaction(Species,scheme)), 
#   #           linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              linetype = 1,color="grey20",size=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(Delta*psi*" (MPa)"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   xlim(-6,0)+
#   ylim(0,1.8)+
#   NULL
# 
# 
# p_dpsi_no_accl <- df_sim_Qi_no %>%
#   rbind(df_sim_Ep_no) %>% 
#   rbind(df_sim_Epo_no) %>%
#   # filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,dpsi,color = scheme),size = 1)+
#   # geom_line(aes(p_leaf,dpsi,color = scheme,group = interaction(Species,scheme)), 
#   #           linetype = 2,size = 1)+
#   # geom_segment(aes(x=psi88S,xend=psi88S, y= 0.5, yend= 0.45,color=scheme),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
#   #              data=df_summary,
#   #              linetype = 1,size=1,
#   #              show.legend = FALSE)+
#   # geom_segment(data=df_summary_sp,aes(x=P50,xend=P50,y=0.5,yend=0.45),
#   #              arrow = arrow(length = unit(0.3, "cm")),
#   #              lineend = c('round'),linejoin = c('mitre'),
# #              linetype = 1,color="grey20",size=1)+
# mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(Delta*psi*" (MPa)"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   xlim(-6,0)+
#   ylim(0,2)+
#   NULL
# 
# ggarrange(p_dpsi_accl,p_dpsi_no_accl, 
#           align='hv', labels=c('a', 'b'),
#           common.legend = T,ncol=1, nrow = 2)
# 
# 
# 
# df_sim_Qi %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   filter(a>0.2) %>% 
#   ggplot()+
#   geom_point(data=df_test,
#              aes(gC,A),
#              color = "grey20",shape=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
#   ylab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   NULL
# 
# 
# 
# 
# 
# 
# df_sim_Qi  %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(var,a/ci,color = scheme),size = 1)+
#   # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#   #           linetype = 2,size = 1)+
#   geom_point(data=df_test,
#              aes(LWP,A/(Ciest*101325/1e6)),
#              color = "grey20",shape=1)+
#   geom_line(aes(gs,a,color = scheme,group = interaction(Species,scheme)),size = 1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("A/"*c[i]~" ("*mu*"mol m"^-2~"s"^-1*"Pa"^-1~")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   NULL
# 
# 
# df_sim_Qi  %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(gs,a/ci,color = scheme),size = 1)+
#   # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#   #           linetype = 2,size = 1)+
#   geom_point(data=df %>% 
#                filter(Species%in%c("Eucalyptus populnea",
#                                    "Eucalyptus pilularis",
#                                    "Quercus ilex")),
#              aes(gC,A/(Ciest*101325/1e6)),
#              color = "grey20",shape=1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(g[s]*" (mol m"^-2~"s"^-1*")"))+
#   ylab(expression(atop("A/"*c[i]~" ("*mu*"mol m"^-2~"s"^-1*"Pa"^-1~")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   NULL
# 
# df_sim_Qi  %>% 
#   rbind(df_sim_Ep) %>% 
#   rbind(df_sim_Epo) %>% 
#   filter(a>0.2) %>%
#   ggplot()+
#   geom_line(aes(a,ci,color = scheme),size = 1)+
#   # geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#   #           linetype = 2,size = 1)+
#   geom_point(data=df %>% 
#                filter(Species%in%c("Eucalyptus populnea",
#                                    "Eucalyptus pilularis",
#                                    "Quercus ilex")),
#              aes(A,(Ciest*101325/1e6)),
#              color = "grey20",shape=1)+
#   mytheme()+
#   labs(color = "")+
#   ylab(expression(c[i]*" (Pa)"))+
#   xlab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1~")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   NULL

################################################################################
# 
# df_sim_Q_low %>%
#   # rbind(df_sim_Ep) %>% 
#   # rbind(df_sim_Epo) %>% 
#   group_by(scheme,Species) %>% 
#   mutate(
#     rprima1 = a/(lead(ci)-lag(ci)),
#     r1 = 1/gs,
#     r2 = 1/lead(gs),
#     rprima2 = lead(rprima1),
#     A1 = a*1e-6,
#     A2 = lead(a)*1e-6,
#     L = (r2-r1)/(A2-A1),
#     sens_g = -(0.5*((A1/(rprima1+r1))+(A2/(rprima2+r2)))*L)) %>% 
#   filter(sens_g<=1,gs>=gs0*0.12,sens_g>=0) %>% 
#   ggplot()+
#   geom_line(aes(var,sens_g,color = scheme,group= scheme),size = 1)+
#   geom_vline(data=df_summary_sp, aes(xintercept=P50),linetype=3)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression("Relative stomatal limitation"))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   guides(colour = guide_legend(nrow = 1))+
#   facet_wrap(~Species)+
#   NULL







# df %>% 
#   ggplot(aes(LWP,gC/A,color=Species))+
#   geom_point(shape = 1)+
#   # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
#   #             se = FALSE,method="lm",linetype = 2, size=0.5)+
#   geom_smooth(se = FALSE, method="gam",formula = y ~ s(x, bs = "cs",k=3))+
#   # facet_wrap(~acclimation)+
#   mytheme2()+
#   theme(
#     legend.title = element_blank(),
#     legend.position = "top",
#     legend.background = element_rect(fill="transparent"))+
#   ylab(expression(atop("WUE"[i]~"("*mu*"mol mol"^{-1}*")")))+
#   xlab(expression(psi*" (MPa)"))+
#   guides(colour = guide_legend(nrow = 2))+
#   NULL




#### EXPLORATORY ALPHA ANALYSIS ####
sla_alpha_mod <- lmerTest::lmer(alpha~SLA..cm2.g.1.+
                                  (1|scheme) + (1|source), 
                                data = df_param_kmax_alpha, weights = log(n_dist)
)
anova(sla_alpha_mod)
sla_alpha <- summary(sla_alpha_mod)$coefficients
sla_sim <- seq(min(df_param_kmax_alpha$SLA..cm2.g.1.,na.rm = TRUE),
               max(df_param_kmax_alpha$SLA..cm2.g.1.,na.rm = TRUE),
               1)
alpha_sim_sla <- sla_alpha[1,1]+sla_alpha[2,1]*sla_sim
alpha_sd_max_sla <- sla_alpha[1,1] + sla_alpha[1,2] + (sla_alpha[2,1]+sla_alpha[2,2])*sla_sim
alpha_sd_min_sla <- sla_alpha[1,1] - sla_alpha[1,2] + (sla_alpha[2,1]-sla_alpha[2,2])*sla_sim
(p1 <- ggplot()+
    stat_summary(data = df_param_kmax_alpha,
                 mapping=aes(SLA..cm2.g.1.,alpha,group=interaction(species,source),size=log(n_dist)),
                 fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="grey20", show.legend = FALSE)+
    geom_line(aes(x = sla_sim, y = alpha_sim_sla),color = "grey20")+
    geom_line(aes(x = sla_sim, y = alpha_sd_max_sla),color = "grey20",linetype=2)+
    geom_line(aes(x = sla_sim, y = alpha_sd_min_sla),color = "grey20",linetype=2)+
    # annotate(geom = "text",x=5,y=0.13,hjust=0,vjust=0,label="***",size=6)+
    mytheme2()+
    scale_radius(range=c(0.2,1))+
    xlab(expression("SLA  [cm"^2~"g"^-1*"]"))+
    ylab(expression(alpha)))

p50_alpha_mod <- lmerTest::lmer(alpha~P50 + 
                                  (1|scheme) + (1|source), 
                                data = df_param_kmax_alpha, weights = log(n_dist)
)
anova(p50_alpha_mod)
p50_alpha <- summary(p50_alpha_mod)$coefficients
p50_sim <- seq(min(df_param_kmax_alpha$P50,na.rm = TRUE),
               max(df_param_kmax_alpha$P50,na.rm = TRUE),
               1)
alpha_sim_p50 <- p50_alpha[1,1]+p50_alpha[2,1]*p50_sim
alpha_sd_max_p50 <- p50_alpha[1,1] + p50_alpha[1,2] + (p50_alpha[2,1]+p50_alpha[2,2])*p50_sim
alpha_sd_min_p50 <- p50_alpha[1,1] - p50_alpha[1,2] + (p50_alpha[2,1]-p50_alpha[2,2])*p50_sim
(p2 <- ggplot()+
    stat_summary(data = df_param_kmax_alpha,
                 mapping=aes(P50,alpha,group=interaction(species,source),size=log(n_dist)),fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="grey20", show.legend = FALSE)+
    geom_line(aes(x = p50_sim, y = alpha_sim_p50),color = "grey20")+
    geom_line(aes(x = p50_sim, y = alpha_sd_max_p50),color = "grey20",linetype=2)+
    geom_line(aes(x = p50_sim, y = alpha_sd_min_p50),color = "grey20",linetype=2)+
    # annotate(geom = "text",x=-7.5,y=0.13,hjust=0,vjust=0,label="**",size=6)+
    mytheme2()+
    scale_radius(range=c(0.2,1))+
    xlab(expression(psi[50]~"[MPa]"))+
    ylab(expression(alpha))
)


# p50_opt_alpha_mod <- lmerTest::lmer(alpha~log(-p50_opt) + 
#                                   (1|scheme) + (1|source/Species), 
#                                 data = df_param_kmax_alpha, weights = log(n_dist)
# )
# anova(p50_opt_alpha_mod)
# p50_opt_alpha <- summary(p50_opt_alpha_mod)$coefficients
# p50_opt_sim <- seq(min(df_param_kmax_alpha$p50_opt,na.rm = TRUE),
#                max(df_param_kmax_alpha$p50_opt,na.rm = TRUE),
#                1)
# alpha_sim_p50_opt <- p50_opt_alpha[1,1]+p50_opt_alpha[2,1]*p50_opt_sim
# alpha_sd_max_p50_opt <- p50_opt_alpha[1,1] + p50_opt_alpha[1,2] + (p50_opt_alpha[2,1]+p50_opt_alpha[2,2])*p50_opt_sim
# alpha_sd_min_p50_opt <- p50_opt_alpha[1,1] - p50_opt_alpha[1,2] + (p50_opt_alpha[2,1]-p50_opt_alpha[2,2])*p50_opt_sim
# (p20 <- ggplot()+
#     geom_point(data = df_param_kmax_alpha,
#                  mapping=aes(p50_opt,alpha,size=log(n_dist)),fun.data=mean_sdl, color="grey20", show.legend = FALSE)+
#     geom_line(aes(x = p50_opt_sim, y = alpha_sim_p50_opt),color = "grey20")+
#     geom_line(aes(x = p50_opt_sim, y = alpha_sd_max_p50_opt),color = "grey20",linetype=2)+
#     geom_line(aes(x = p50_opt_sim, y = alpha_sd_min_p50_opt),color = "grey20",linetype=2)+
#     # annotate(geom = "text",x=-7.5,y=0.13,hjust=0,vjust=0,label="**",size=6)+
#     mytheme2()+
#     scale_radius(range=c(0.2,1))+
#     xlab(expression(psi[50]~"[MPa]"))+
#     ylab(expression(alpha))
# )


hv_alpha_mod <- lmerTest::lmer(alpha~log(Huber.value)+  
                                 (1|scheme) + (1|source), 
                               data = df_param_kmax_alpha, weights = log(n_dist))
anova(hv_alpha_mod)
hv_alpha <- summary(hv_alpha_mod)$coefficients
hv_sim <- seq(min(log(df_param_kmax_alpha$Huber.value),na.rm = TRUE),
              max(log(df_param_kmax_alpha$Huber.value),na.rm = TRUE),
              0.1)
alpha_sim_hv<- hv_alpha[1,1]+hv_alpha[2,1]*hv_sim
alpha_sd_max_hv <- hv_alpha[1,1] + hv_alpha[1,2] + (hv_alpha[2,1]+hv_alpha[2,2])*hv_sim
alpha_sd_min_hv <- hv_alpha[1,1] - hv_alpha[1,2] + (hv_alpha[2,1]-hv_alpha[2,2])*hv_sim
(p3 <- ggplot()+
    stat_summary(data = df_param_kmax_alpha,
                 mapping=aes(log(Huber.value),alpha,group=interaction(Species,source),size=log(n_dist)),
                 fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="grey20", show.legend = FALSE)+
    geom_line(aes(x = hv_sim, y = alpha_sim_hv),color = "grey20")+
    geom_line(aes(x = hv_sim, y = alpha_sd_max_hv),color = "grey20",linetype=2)+
    geom_line(aes(x = hv_sim, y = alpha_sd_min_hv),color = "grey20",linetype=2)+
    annotate(geom = "text",x=-7,y=0.13,hjust=0,vjust=0,label="**",size=6)+
    mytheme2()+
    scale_radius(range=c(0.2,1))+
    xlab(expression("Huber value (ln("*cm[sw]^2~m[leaf]^-2*"))"))+
    ylab(expression(alpha))
)


k_alpha_mod <- lmerTest::lmer(alpha~log(K.scale)+
                                (1|scheme) + (1|source/Species), 
                              data = df_param_kmax_alpha, weights = log(n_dist)
)
anova(k_alpha_mod) 
k_alpha <- summary(k_alpha_mod)$coefficients
k_sim <- seq(min(log(df_param_kmax_alpha$K.scale),na.rm = TRUE),
             max(log(df_param_kmax_alpha$K.scale),na.rm = TRUE),
             0.1)
alpha_sim_k <- k_alpha[1,1]+k_alpha[2,1]*k_sim
alpha_sd_max_k <- k_alpha[1,1] + k_alpha[1,2] + (k_alpha[2,1]+k_alpha[2,2])*k_sim
alpha_sd_min_k <- k_alpha[1,1] - k_alpha[1,2] + (k_alpha[2,1]-k_alpha[2,2])*k_sim
(p4 <- ggplot()+
    stat_summary(data = df_param_kmax_alpha,
                 mapping=aes(log(K.scale),alpha,group=interaction(species,source),size=log(n_dist)),
                 fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="grey20", show.legend = FALSE)+
    # stat_summary(data = df_param_kmax_alpha,
    #              aes(log(55.55*KL..kg.m.1.MPa.1.s.1.),alpha,
    #                  group=interaction(species,source)),
    #              fun.data=mean_sdl, fun.args = list(mult=1),
    #              geom="pointrange", shape=1)+
    geom_line(aes(x = k_sim, y = alpha_sim_k),color = "grey20")+
    geom_line(aes(x = k_sim, y = alpha_sd_max_k),color = "grey20",linetype=2)+
    geom_line(aes(x = k_sim, y = alpha_sd_min_k),color = "grey20",linetype=2)+
    # annotate(geom = "text",x=-3.1,y=0.13,hjust=0,vjust=0,label="*",size=6)+
    mytheme2()+
    scale_radius(range=c(0.2,1))+
    xlab(expression(K[leaf]*" [ln("*mol[H2O]~m^-1~MPa^-1~s^-1*")]"))+
    ylab(expression(alpha))
)
b_alpha_mod <- lmerTest::lmer(alpha~log(b) + #(1|Species)+ 
                                (1|scheme) + (1|source), 
                              data = df_param_kmax_alpha, weights = log(n_dist))
anova(b_alpha_mod)
b_alpha <- summary(b_alpha_mod)$coefficients
b_sim <- seq(min(log(df_param_kmax_alpha$b),na.rm = TRUE),
             max(log(df_param_kmax_alpha$b),na.rm = TRUE),
             1)
alpha_sim_b <- b_alpha[1,1]+b_alpha[2,1]*b_sim
alpha_sd_max_b <- b_alpha[1,1] + b_alpha[1,2] + (b_alpha[2,1]+b_alpha[2,2])*b_sim
alpha_sd_min_b <- b_alpha[1,1] - b_alpha[1,2] + (b_alpha[2,1]-b_alpha[2,2])*b_sim
(p5 <- ggplot()+
    stat_summary(data = df_param_kmax_alpha,
                 mapping=aes(log(b),alpha,group=interaction(Species,source),size=log(n_dist)),
                 fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="grey20", show.legend = FALSE)+
    # geom_line(aes(x = b_sim, y = alpha_sim_b),color = "grey20")+
    # geom_line(aes(x = b_sim, y = alpha_sd_max_b),color = "grey20",linetype=2)+
    # geom_line(aes(x = b_sim, y = alpha_sd_min_b),color = "grey20",linetype=2)+
    # annotate(geom = "text",x=0.1,y=0.13,hjust=0,vjust=0,label="***",size=6)+
    mytheme2()+
    scale_radius(range=c(0.2,1))+
    xlab(expression("ln(b)"))+
    ylab(expression(alpha))
)

Iabs_alpha_mod <- lmerTest::lmer(alpha~ Iabs_growth +# (1|)+ 
                                   (1|scheme) + (1|source/Species), 
                                 data = df_param_kmax_alpha, weights = log(n_dist))
anova(Iabs_alpha_mod)
Iabs_alpha <- summary(Iabs_alpha_mod)$coefficients
Iabs_sim <- seq(min(df_param_kmax_alpha$Iabs_growth,na.rm = TRUE),
                max(df_param_kmax_alpha$Iabs_growth,na.rm = TRUE),
                1)
alpha_sim_Iabs <- Iabs_alpha[1,1]+Iabs_alpha[2,1]*Iabs_sim
alpha_sd_max_Iabs <- Iabs_alpha[1,1] + Iabs_alpha[1,2] + (Iabs_alpha[2,1]+Iabs_alpha[2,2])*Iabs_sim
alpha_sd_min_Iabs <- Iabs_alpha[1,1] - Iabs_alpha[1,2] + (Iabs_alpha[2,1]-Iabs_alpha[2,2])*Iabs_sim
(p6 <- ggplot()+
    geom_point(data = df_param_kmax_alpha,
               mapping=aes(Iabs_growth,alpha,group=interaction(Species,source),size=log(n_dist)), 
               color="grey20", show.legend = FALSE)+
    # geom_line(aes(x = Iabs_sim, y = alpha_sim_Iabs),color = "grey20")+
    # geom_line(aes(x = Iabs_sim, y = alpha_sd_max_Iabs),color = "grey20",linetype=2)+
    # geom_line(aes(x = Iabs_sim, y = alpha_sd_min_Iabs),color = "grey20",linetype=2)+
    # annotate(geom = "text",x=0.1,y=0.13,hjust=0,vjust=0,label="***",size=6)+
    mytheme2()+
    scale_radius(range=c(0.8,4))+
    xlab(expression(I[growth]~"["*mu*"mol"~m^-2~s^-1*"]"))+
    ylab(expression(alpha))
)


ggarrange(p1,p2,p3,p4,p5,p6,
          align='hv', labels=c('a', 'b','c','d','e','f'),
          ncol=2, nrow = 3)



################################################################################


# vcmax_jmax_params <- df_a_fix_param%>% 
#   left_join(df_a_fix %>% 
#               dplyr::select(Species,ca,T,D,Iabs_used) %>% 
#               group_by(Species) %>% 
#               summarise_all(mean, na.rm=TRUE)) %>% 
#   filter(Species%in%c("Picea abies","Populus tremula"),
#          acclimation == "Acclimated",
#          dpsi == FALSE)
# 
# 
# df_sim_jmax <- vcmax_jmax_params %>%
#   rowwise() %>%
#   do(mapped = get_acclimated_response(.,tc_now=.$T,ppfd_now=.$Iabs_used,vpd_now=.$D*101325,co2_now = .$ca)) %>% 
#   unnest(cols = c(mapped))
# df_summary_jmax <- df_sim_jmax  %>% 
#   group_by(Species,scheme) %>% 
#   dplyr::select(psi88S) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# df_test_jmax <- df_a_fix %>% 
#   filter(Species%in%c("Picea abies","Populus tremula"))
# 
# df_summary_jmax <- df_test_jmax %>% 
#   group_by(Species) %>% 
#   dplyr::select(P50) %>% 
#   summarise_all(mean,na.rm=TRUE)
# 
# 
# df_summary_jmax_error <- df_a_fix %>% filter(Species%in%c("Picea abies","Populus tremula"),
#               acclimation == "Acclimated") %>% 
#   dplyr::select(LWP,Species,week,jmax_obs) %>% 
#   group_by(Species,week) %>% 
#   drop_na() %>% 
#   summarise(
#     jmax_sd = sd(jmax_obs,na.rm=TRUE),
#     LWP_sd = sd(LWP,na.rm=TRUE),
#     jmax_obs=mean(jmax_obs,na.rm=TRUE),
#     LWP=mean(LWP,na.rm=TRUE)) %>% 
#   drop_na()
# 
# df_sim_jmax %>%
#   ggplot()+
#   geom_line(aes(var,jmax,color = scheme),size = 1)+
#   geom_line(aes(p_leaf,jmax,color = scheme,group = interaction(Species,scheme)), 
#             linetype = 2,size = 1)+
#   geom_point(data = df %>% filter(Species%in%c("Picea abies","Populus tremula"),
#                                   acclimation == "Acclimated") %>%
#                group_by(Species,LWP) %>% summarise(jmax_obs=mean(jmax_obs,na.rm=TRUE)),
#              mapping=aes(LWP,jmax_obs),
#              color="grey50")+
#   geom_point(data=df_summary_jmax_error,mapping = aes(LWP,jmax_obs))+
#   geom_errorbarh(aes(y=jmax_obs,
#                      xmin=LWP-LWP_sd,
#                      xmax=LWP+LWP_sd),
#                  data=df_summary_jmax_error,
#                  height=2)+
#   geom_errorbar(aes(x=LWP,
#                     ymin=jmax_obs-jmax_sd,
#                     ymax=jmax_obs+jmax_sd),
#                 data=df_summary_jmax_error,
#                 width=0.1)+
#   mytheme()+
#   labs(color = "")+
#   xlab(expression(psi*" (MPa)"))+
#   ylab(expression(atop("J"[max]*" ("*mu*"mol m"^-2~"s"^-1*")")))+
#   scale_colour_manual(breaks = col_df$scheme, 
#                       values = unique(as.character(col_df$col)))+
#   theme(legend.position = "top",
#         legend.background = element_rect(fill="transparent"))+
#   facet_wrap(~Species)+
#   xlim(-7,0)+
#   guides(colour = guide_legend(nrow = 1))+
#   NULL


