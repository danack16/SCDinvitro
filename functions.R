

del.V2<- function(df){
  ## deletes V2 column
  df$V2<- NULL
  df
}

# creating function that interprets lipid name and breaks it down
# uses different strings for different lipids
CL.split<- function(df, cl){
  if(cl=="CL1"){
    chains<- cl1_ch
    nms<- c("Label", "Species", "ch1", "ch1DBs")
  } 
  if(cl=="CL2"){
    chains<- cl2_ch
    nms<- c("Label", "Species", "Add1", "ch1", "ch1DBs","Add2", "ch2", "ch2DBs")
  }
  if(cl=="CL3"){
    chains<- cl3_ch
    nms<- c("Label", "Species", "ch1", "ch1DBs", "ch2", "ch2DBs", "ch3", "ch3DBs")
    ##why do CL1 and CL3 not contain Add1 column?
  }
  #else stop("unknown chain number") # doesn't work with this 
  
  FA_chains<- as.data.frame(str_match(df$Label, chains))
  colnames(FA_chains)<- nms
  merge(df, FA_chains, by="Label")
}

tot.unsat2<- function(df, cl){
  ## groups lipids into bins based on number of double bonds
  ## calculates sum of each bin to get total abundance of lipids with 0 DBs, etc
  if (!"ch1DBs" %in% names(df)){
    stop("df must contain ch1DBs column")
  }
  exclude<- c('Label',"ch1", "ch1DBs", "Species")
  if(cl=="CL1"){
    df$unsat<- as.numeric(as.character(df$ch1DBs))
  }
  if(cl=="CL2"){
    df$unsat<- as.numeric(as.character(df$ch1DBs))+as.numeric(as.character(df$ch2DBs))
    exclude<- c(exclude, c("ch2", "ch2DBs", "Add1", "Add2"))
  }
  if(cl=="CL3"){
    df$unsat<- as.numeric(as.character(df$ch1DBs))+as.numeric(as.character(df$ch2DBs))+ as.numeric(df$ch3DBs)
    exclude<- c(exclude, c("ch2", "ch2DBs", "ch3", "ch3DBs"))
  }
  df2<-setDT(df)[,lapply(.SD, sum), by=unsat, .SDcols=-exclude]
  df2
}

scal<- function(DAT){
  ##scales each sample, so that total abundance = 1
  ## this allows for the evaluation of changes in the fraction of lipids that are 
  ## relativels saturated in downstream steps
  x<-DAT$unsat
  nms<- colnames(DAT)
  df2<-sweep(DAT,2,colSums(DAT), "/")
  df2$unsat<- x
  colnames(df2)<- nms
  df2[order(df2$unsat),]
}

scal2<- function(df, exclude = c()){
  ##df must contain column named "unsat"
  ## same as scal() function, except that it allows for exclusion of columns
  ## this makes it easier to use this function on dataframes in which some columns 
  ## contain information other than abundance
  exclude2<- c("unsat", exclude)
  excluded<- select(df, one_of(exclude2))
  included<- names(df)[!(names(df)%in% exclude2)]
  DAT<- select(df, one_of(included))
  df2<-sweep(DAT,2,colSums(DAT), "/")
  df3<- cbind(excluded, df2)
  df3[order(df3$unsat),]
}

centr<- function(df, controls){
  # takes a dataframe and a vector of sample names
  # calculates mean of each saturatio value across control samples
  # subtracts mean from all values
  df_ctrl<- df[,names(df) %in% controls]
  df2<- sweep(df, 1, rowMeans(df_ctrl), "-")
  df2$unsat<- df$unsat
  df2
}

diff.plot<- function(df, exclude=c(), order=c()){
  ## function that takes centered dataframe created by centr()
  ## converts into long format
  ## graphs it in ggplot2
  ## contains argument that takes list of column names to exclude
  ## contains argument that allows user to decide which order to graph in (levels for ggplot2)
  ## naming of samples has to be as follows:
  ### column names match group names with "_1", "_2", etc to identify individual samples
  
  long_df<- gather(df, Sample, Measurement, 2:length(colnames(df)))
  long_df$Group<-unlist(strsplit(as.character(long_df$Sample), '[1-5]$'))
  long_df$Group<- gsub("_$", "", long_df$Group)
  long_df2<- long_df[!(long_df$Sample %in% exclude),]
  try(long_df2$Group<- factor(long_df2$Group, levels = order))
  ##plotting  
  p<-ggplot(long_df2, 
            aes(x=factor(unsat), y=Measurement, fill= factor(Group)))
  p<- p +    geom_boxplot()
  p<- p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  labs(x="Double Bonds")
  p
}

mneffect<- function(df, control, treatment){
  #function that graphs changes in saturation profile between two groups
  #takes a dataframe and two strings of column names
  #output is gg barplot
  
  #### include some stop options if control and treatment columns are not present in names(df)
  
  df_control<- df[, names(df) %in% control]
  df_treat<- df[, names(df) %in% treatment]
  df_diff<- sweep(df_treat, 1, rowMeans(df_control), "-")
  
  ##Create df with all the data required for ggplot
  df_diff2<- data.frame(DBs=df[,1])
  df_diff2$mean<- apply(df_diff,1,mean)
  df_diff2$sd<- apply(df_diff,1,se)
  df_diff2$ctrlmn<- rowMeans(df_control)
  
  
  
  z<- ggplot(df_diff2, aes(x=DBs, y=mean, fill=mean>0))+
    geom_bar(stat="identity", position="identity", color="black", size=0.5)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_x_continuous(breaks=seq(0,9,1))+
    scale_y_continuous(breaks=seq(-0.1,0.3,0.05))+
    scale_fill_manual(values=c("#CCEEFF", "#FFDDDD"), guide=FALSE) +
    geom_hline(aes(yintercept=0))+labs(y="relative change")+
    coord_cartesian(xlim = c(0.5, 8.5), ylim=c(-.2, 0.45))+
    theme(axis.title.x=element_text(size = 25), axis.text.x=element_text(size=30),
          axis.title.y=element_text(size = 25), axis.text.y=element_text(size=30))+
    theme_bw()+
    
    geom_point(aes(size=ctrlmn, x= DBs, y=0.35), 
               shape = 21, color= "black", fill= "cornsilk")+
    scale_size_area(max_size=6.5, guide=FALSE)+
    annotate("text", x=4.5, y=0.43, label= "TG species", size=5)+
    annotate("text", x=4.5, y=0.39, label= "Abundance", size=3.5)
  z  
  
}

mneffect2<- function(df, control, treatment){
  #function that graphs changes in saturation profile between two groups
  #takes a dataframe and two strings of column names
  #output is gg barplot
  
  #### include some stop options if control and treatment columns are not present in names(df)
  
  df_control<- df[, names(df) %in% control]
  df_treat<- df[, names(df) %in% treatment]
  df_diff<- sweep(df_treat, 1, rowMeans(df_control), "-")
  
  ##Create df with all the data required for ggplot
  df_diff2<- data.frame(DBs=df[,1])
  df_diff2$mean<- apply(df_diff,1,mean)
  df_diff2$sd<- apply(df_diff,1,se)
  df_diff2$ctrlmn<- rowMeans(df_control)
  
  
  
  z<- ggplot(df_diff2, aes(x=DBs, y=mean, fill=mean>0))+
    geom_bar(stat="identity", position="identity", color="black", size=0.5)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
  z  
  
}



TG_effects<- function(lst, grp1, grp2){
  ## takes two text strings, identifies columns in df with those names
  ## takes list, extracts "TG" element as dataframe
  ## performs mneffect() function on those columns
  
  df<-lst[["TG"]]
  cols1<- grep(grp1, names(df), value=TRUE)
  cols2<- grep(grp2, names(df), value=TRUE)
  p<- mneffect2(df, cols1, cols2)
  p
}

FAsplit<- function(df){
  ## function to take dataframe of lipids and split into individual FA chains
  ## df needs to contain columns with chain length and number of double bonds for each FA
  ## this was done above using stringr breakup of lipid name into constituents
  ## takes each TG and introduces that entry three times, once for each FA chain, 
  ## preserving only chain length/saturation for one of the chains in each entry. 
  ## assumes the following naming for columns: 
  ### All sample column names include "_"
  ### chain length columns named "ch1", "ch2", "ch3" ("ch2" and "ch3" might be absent)
  ### DB number columns "ch1DBs", "ch2DBs", "ch3DBs" (last two might be absent)
  
  ch1<- c(grep("ch1$", names(df)), grep("ch1DBs", names(df)))
  ch2<- c(grep("ch2$", names(df)), grep("ch2DBs", names(df)))
  ch3<- c(grep("ch3$", names(df)), grep("ch3DBs", names(df)))
  
  if (!("ch1" %in% names(df) & "ch1DBs" %in% names(df))){
    stop("chain length and saturation not provided in correct format")
  }
  cols1<- c(ch1, grep("_", names(df)))
  df2<- df[,cols1]
  names(df2)<- c("length", "unsat", names(df)[grep("_", names(df))]) 
  
  if ("ch2" %in% names(df) & "ch2DBs" %in% names(df)){
    cols2<- c(ch2, grep("_", names(df)))
    df_ch2<- df[,cols2]
    names(df_ch2)<- c("length", "unsat", names(df)[grep("_", names(df))]) 
    df2<- rbind(df2, df_ch2)
  }
  
  if ("ch3" %in% names(df) & "ch3DBs" %in% names(df)){
    cols3<- c(ch3, grep("_", names(df)))
    df_ch3<- df[,cols3]
    names(df_ch3)<- c("length", "unsat", names(df)[grep("_", names(df))]) 
    df2<- rbind(df2, df_ch3)
  }
  df2
}

desat.index<- function(df, exclude=c()){
  ## function takes df with following columns
  ### Samples with lipid quantity measurements
  ### each row is one lipid measurement
  ### column named "unsat" containing number of double bonds in that lipid
  ## calculates desaturation index for each sample
  ## allows for exclusion of columns from calculation
  
  exclude2<- c("unsat", exclude)  #vector of columns to exclude from calculation
  excluded<- select(df, one_of(exclude2))   # df of excluded columns
  include<- names(df)[!(names(df)%in% exclude2)]  # vector of columns to include
  DAT<- select(df, one_of(include))  #df of columns of interest
  excluded$unsat<- as.numeric(excluded$unsat)
  DAT2<- sweep(DAT, 1, excluded$unsat, "*")
  Vec1<- colMeans(DAT2)
  names(Vec1)<- names(DAT2)
  Vec1
}

desat.effect<- function(df, control, treatment, lipids){
  # calculates change in desaturation index between a control and a treatment group
  # plots boxplot for selected lipids 
  df_control<- filter(df, Sample %in% control) %>% select(one_of(lipids))
  df_treat<- filter(DesatIndex2, Sample %in% treatment) %>% select(one_of(lipids))
  df_diff<- sweep(df_treat, 2, colMeans(df_control), "-")
  long<- gather(df_diff)
  long$value<- long$value*1000
  long<- filter(long, key %in% lipids)
  p<- ggplot(long, aes(x=key, y=value))+ geom_boxplot()
  p
}

mk.names<- function(string, df){
  # uses grep string to identify all column names in a dataframe that contain each string in vector "string"
  # takes dataframe and vector of character strings corresponding to groupings
  # output is a list. name of list elements corresponds to group, contents to samples in that group
  
  out1<- vector("list", length = length(string))
  for (i in 1:length(string)){
    x <-  names(df)[grep(string[i], names(df))]
    out1[[i]]<- x
    names(out1)[i]<- string[i]
  }
  out1
  
}

abund<- function(df, exclude=c()){
  ## df: dataframe containg raw abundance for each lipid (each row = 1 lipid, each column= one sample). 
  ## exclude: columns to exclude from output. Auto excludes "Label" column, which contains lipid names
  ## function operation
  ### take column means and copy into a dataframe
  #### structure of this dataframe: 
  #### columns are replicate names
  #### each row corresponds to different lipid type
  ###output this dataframe
  exclude2<- c("Label", exclude)
  excluded<- select(df, one_of(exclude2))
  included<- names(df)[!(names(df)%in% exclude2)] #vector of included columns
  df2<- select(df, one_of(included))
  colSums(df2)
  # df2 
}

group.mean<- function(vec, df){
  ## takes dataframe, subsets only columns included in "vec"
  ## calculates rowMeans for reduced dataframe
  df2<- df[, colnames(df) %in% vec]
  mns<- rowMeans(df2, na.rm = TRUE)
  mns
}

ht.map<- function(names, df, ctrl_group){
  ## df containing each sample in a differnet row, with total abundance of each lipid type in each column
  ## names: list of vectors corresponding to sample groups. name of vector is group name, elements of vector are replicates of that group. 
  ### each dataframe, for loop over elements of "names"
  ### identify columns corresponding to elements of names vector
  
  df2<- vapply(names, group.mean, df, FUN.VALUE = numeric(nrow(df)))
  if (!ctrl_group %in% colnames(df2)){
    stop("Incorrect group name for control group")
  }
  ctrl<- df2[,colnames(df2)==ctrl_group]
  df3<- sweep(df2, 1, ctrl, "/")
  df4<- as.data.frame(df3) #as.data.frame(abundance4[, 2:ncol(abundance4)])
  df4$Lipid<- c(rownames(df4))
  long_abundance<- gather(df4, key = Group, value= Measurement, 1:(ncol(df4)-1))
  long_abundance
}

samp.cens<- function(df, to_censor){
  ## take df with samples in columns according to sample name and each row = different lipid
  ## take vector of sample names that need to be removed
  ## output dataframe with values in those columns replaced by NA
  if (!(all(to_censor %in% names(df)))){
    stop("Incorrect sample names")
  }
  df[,to_censor]<- NA
  df
  #cens<- setdiff(names(df), to_censor)
  #   for (i in 1:length(to_censor)){
  #     
  #   }
}

sem<- function(vec){
  sm<- sd(vec)/sqrt(length(vec))
}


lip.abund<- function(lip,df, groups, nms= names1){
  ## takes dataframe, subsets groups = vector of character strings corresponding to groups of interest within names. 
  ### dataframe: rows= total abundance by lipid type, columns = samples
  ## takes names, a list of character strings. Element names = group names, element contents = samples within that group
  ## takes lip, a character string corresponding to one of the following lipids present as row in df
  ## extracts abundance of lipid selected through "lip" for groups selected in "groups"
  ## outputs ggplot object
  
  samples<- c()
  names2<- nms[groups]
  for(i in 1:length(groups)){
    samples<- c(samples, names2[[i]])
  }
  df2<- select(df, one_of(samples))%>% add_rownames(var = "lipid")%>%filter(lipid==lip)%>% gather(Sample, Measurement, 2:(length(samples)+1))
  df2$Group<- NA
  for (i in 1:length(groups)){
    #     grp<-groups[i]
    #     gr_samples<- 
    df2$Group[df2$Sample %in% names2[[i]]]<- names(names2)[i]
  }
  
  df3<- df2%>%group_by(Group)%>% mutate(lQntl = quantile(Measurement, probs = 0.25, na.rm = TRUE),
                                        uQntl = quantile(Measurement, probs = 0.75, na.rm = TRUE), lBound = lQntl - 1.5 * (uQntl - lQntl),
                                        uBound = uQntl + 1.5 * (uQntl - lQntl))
  df3$OL<-df3$Sample
  df3$OL[!(df3$Measurement<=df3$lBound |df3$Measurement>=df3$uBound)]<- NA
  
  
  ggplot(df3, aes(x=Group, y=Measurement))  +geom_boxplot() +geom_text(aes(label = OL), size = 3,
                                                                       position = position_jitter(width = 0.1))
}

se <- function(x) sqrt(var(x)/length(x))