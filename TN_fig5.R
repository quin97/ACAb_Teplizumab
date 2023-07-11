# https://github.com/rdboyes/forester
library(forester) 
df2<-data.frame(Subgroup = c("Treatment effect in entire cohort","Baseline B.longum - IgG2", "  Low", "  High",
                             "Baseline E.faecalis - IgG2", "  Low", "  High",
                             "DR4", "  DR4-Absent", "    Baseline D.invisus - IgG2 Low","    Baseline D.invisus - IgG2 High", "  DR4-Present",
                             "    Baseline D.invisus - IgG2 Low","    Baseline D.invisus - IgG2 High"),
                No.of.Participant = c(62,NA, 39, 23, NA, 29, 33, NA, 21, 11,10, 41,20,21),
                Estimate = c(-0.68878, NA, -0.19115, -1.71636, NA, -0.24962, -0.85320, NA,0.77087,-0.84001,1.2564,-1.62999,-1.5832,-2.90287),
                CI.low = c(-1.03611,NA, -0.74988, -2.32079, NA,-0.91563,-1.38659,NA,0.07614,-2.11709,-0.5633,-2.10912,-2.8113,-3.97423),
                CI.high = c(-0.34145,NA, 0.36758, -1.11193, NA, 0.41639,-0.31981,NA,1.4656,0.43707,3.0761,-1.15086,-0.3551,-1.83151))
colnames(df2)<-c("Subgroup","No. of Participant","Estimate log(HR)","CI low","CI high")

forester(left_side_data = df2[,1:2],
         estimate = df2$Estimate,
         ci_low = df2$`CI low`,
         ci_high = df2$`CI high`,
         nudge_y = .2,
         display = F,
         xlim = c(-4, 3.0),
         null_line_at = 0,
         arrows = TRUE, 
         arrow_labels = c("Teplizumab Better", "Placebo Better"),
         estimate_precision = 2,
         estimate_col_name = "log(Hazard Ratio) (95% CI)",
         justify = c(0,0.5,0,1),
         render_as = "svg",
         file_path = "forest2log.svg")