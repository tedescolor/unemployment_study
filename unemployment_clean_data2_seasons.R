if (!require("pacman")) install.packages("pacman")
pacman::p_load("quantreg","foreach","Rlab","optimbase","survminer","EnvStats",
               "SimDesign","zeallot", "survival","utils","MVN","dppmix",
               "microbenchmark","reshape2","dplyr","gdata","numbers","pracma","doParallel","foreach",
               "abind","icesTAF")


#load data
df.original = read.csv(file = '/home/lorenzo/R/phd/data_master_proef.csv',sep = ';')
set.seed(123)
size.subset = nrow(df.original)
print(paste("percentage used data: ", round(size.subset / nrow(df.original) * 100,1), "%",sep = ""))
df = df.original[sample.int(nrow(df.original),size.subset),]
#df = df.original


#names of the variables
names <- c("region","education","cohort","sex","age",
               "driving","trasportation",
               "status0","status1","status2","status3","status4","status5","status6","status7","status8","status9","status10","status11","status12","status13","status14","status15","status16","status17","status18","status19","status20","status21","status22","status23","status24","status25","status26","status27","status28","status29","status30","status31","status32","status33","status34","status35","status36","status37","status38","status39","status40","status41","status42","status43","status44","status45","status46","status47","status48","status49","status50","status51","status52","status53","status54","status55","status56","status57","status58","status59","status60","status61","status62","status63","status64","status65","status66","status67","status68","status69","status70","status71","status72","status73","status74","status75","status76","status77","status78","status79","status80","status81","status82","status83","status84","status85","status86","status87","status88","status89","status90","status91","status92","status93","status94","status95","status96","status97","status98","status99","status100","status101","status102","status103","status104","status105","status106","status107","status108","status109","status110","status111","status112","status113","status114","status115","status116","status117","status118","status119","status120","status121","status122","status123","status124","status125","status126","status127","status128","status129","status130","status131","status132","status133","status134","status135","status136","status137","status138","status139","status140","status141","status142","status143","status144","status145","status146","status147","status148","status149","status150","status151",
               "knowledge_dutch","immigrant","id")

df = df[!duplicated(df$id),]
names(df) <- names
df$year = as.integer(sapply(df$cohort, substr, start = 1, stop = 4)) # extract year from cohort and convert it into integer
df$month = as.integer(sapply(df$cohort, substr, start = 5, stop = 6)) # same for month
#set the years of comparison
#compare the seasons of unemployment before cut.year(<Jen) and after(>=Jen)
cut.month = 1
#output setting
current.comparison.string = paste(cut.year,"unemployment_study",size.subset,sep = "")  #string to append to every image name before saving
file.output = paste(current.comparison.string,"_output.txt",sep = "") #name file output
# mkdir(current.comparison.string)#create a folder for the output
# setwd(paste("./",current.comparison.string,sep = "")) #enter in that folder


min.year = min(df$year) #min year of inflation
min.month = min(df[df$year == min.year, "month"]) #min month of inflation in min.year
max.year = max(df$year) #same but max 
max.month = max(df[df$year == max.year, "month"]) #same but max
df$inflation =  (df$year -  min.year) * 12 + (df$month -min.month) #count the month of the inflation  counting from month 0 = inflation.year/inflation.month (2007/01)

last.index.season1 = (cut.year -  min.year) * 12 + (cut.month - min.month) - 1
df$last.index.season1  = as.integer(last.index.season1 - df$inflation) #last index of follow up realtive to season 1

#save only the status variables in df.status
df.status = df%>% select(last.index.season1,status0, status1, status2, status3, status4, status5, status6, status7, status8, status9, status10, status11, status12, status13, status14, status15, status16, status17, status18, status19, status20, status21, status22, status23, status24, status25, status26, status27, status28, status29, status30, status31, status32, status33, status34, status35, status36, status37, status38, status39, status40, status41, status42, status43, status44, status45, status46, status47, status48, status49, status50, status51, status52, status53, status54, status55, status56, status57, status58, status59, status60, status61, status62, status63, status64, status65, status66, status67, status68, status69, status70, status71, status72, status73, status74, status75, status76, status77, status78, status79, status80, status81, status82, status83, status84, status85, status86, status87, status88, status89, status90, status91, status92, status93, status94, status95, status96, status97, status98, status99, status100, status101, status102, status103, status104, status105, status106, status107, status108, status109, status110, status111, status112, status113, status114, status115, status116, status117, status118, status119, status120, status121, status122, status123, status124, status125, status126, status127, status128, status129, status130, status131, status132, status133, status134, status135, status136, status137, status138, status139, status140, status141, status142, status143, status144, status145, status146, status147, status148, status149, status150, status151)


# save only the NOT-status variables in df.nostatus
df.nostatus = df %>% select(id,region,education,sex,age,driving,trasportation,knowledge_dutch,immigrant)
head(df.nostatus)


index.unemployed.season1 = apply(df.status ,1, function(row){
  last.index.season1 = row["last.index.season1"]
  status = row[-c(1)] 
  return(any(status[1:last.index.season1] == 1) )
}) # index of rows that presents unemployment status (status = 1) in year1 

index.unemployed.season2 = apply(df.status ,1, function(row){
  first.index.season2 = row["last.index.season1"] + 1
  status = row[-c(1)] 
  status = status[first.index.season2:152]
  become.employed.index = match(2,status)
  if(!is.na(become.employed.index)){ #if there is a significant value
    return(any(status[become.employed.index:length(status)] == 1) ) #check if it RETURN unemployed during season 2
  }
  return(FALSE)

}) 

#take only ineteresting indices
interesting.index = index.unemployed.season1 & index.unemployed.season2 #index that have unemployment status in year1, then find job, than return unemployed in year 2
df.nostatus = df.nostatus[interesting.index,] 
df.status = df.status[interesting.index,]


#for every index find compute the first time/censored variable and the second

dmpstatus.copy = df.status# make a copy of df.status for computatation of the first spell that last from start year1 to start year2
dmpstatus.copy[dmpstatus.copy==2]<-0 #substitute the value 2 (= Working at the end of the month) as 0 (= Status could not be determined)
dmpstatus.copy[dmpstatus.copy==3]<-0 #substitute the value  (= Not in search of a job and not work-ing at the end of the month) as 0 (= Status could not be determined)
 #create df1 as copy of status

df.nostatus$utime1 = apply(dmpstatus.copy, MARGIN = 1, FUN = function(row){
  last.index.season1 = row["last.index.season1"]
  status = row[-c(1)] 
  return(match(0,status[1:last.index.season1])) 
}) - 1 # find the index of first zero (i.e. number of months of eunemployment +1) -1
df.nostatus[,"censored1"] = is.na(df.nostatus$utime1) # NA means still unemployed at the end of the study so censored.
df.nostatus$utime1[is.na(df.nostatus$utime1)] <- df.status$last.index.season1[is.na(df.nostatus$utime1)]# substitute NA value with the censored value (last month of follow up in season 1

# for every row understand if the first zero was for a censoring reason (0, or 3 ) or for employment reason (2). In the first case put censored equal to true
for(j in 1:nrow(df.nostatus)){
  df.nostatus[j,"censored1"] =  df.nostatus[j,"censored1"] |
    df.status[j,df.nostatus[j,"utime1"]+1 +1 ] == 3 | #+1 becouse index inr start from 1 and our index of months starts from 0, +1 becouse we add the inflation first column
    df.status[j,df.nostatus[j,"utime1"]+1 +1 ] == 0
}
table(df.nostatus$censored1)



df.nostatus$utime2 = NA
df.nostatus$censored2 = NA
for(j in 1:nrow(df.nostatus)){
  first.index.season2 = df.status[j,"last.index.season1"] + 1
  status = df.status[j,-c(1)] 
  #print(status)
  status = status[first.index.season2:152] 
  #print(status)
  become.employed.index = match(2,status)
  #print(become.employed.index)
  status = status[become.employed.index:length(status)]
  #print(status)
  return.unemployed.index = match(1,status)  #index of first 0 after first 1 + index of first 1
  #print(return.unemployed.index)
  status = status[return.unemployed.index:length(status)]
  #print(status)
  eventual.employ.again.time = match(c(0,2,3),status)
  #print(eventual.employ.again.time)
  eventual.employ.again.time = eventual.employ.again.time[!is.na(eventual.employ.again.time)] #remove NA
  if(length(eventual.employ.again.time)> 0){#become employed or censored before the end of the study
    #print(min(eventual.employ.again.time))
    again.no.unemployed = min(eventual.employ.again.time)
    df.nostatus[j,"utime2"] = again.no.unemployed -1
    df.nostatus[j,"censored2"] = status[again.no.unemployed] != 2  
    }
  else{
    df.nostatus[j,"utime2"] = 152 -  (first.index.season2 + become.employed.index + return.unemployed.index)
    df.nostatus[j,"censored2"] = TRUE
  } #check if it is AGAIN unemployed during season 2
}
# ##for every index find compute the SECOND time/censored variable and the second
# df.nostatus$utime2 = apply(df.status, MARGIN = 1, FUN = function(row){
#   first.index.season2 = row["last.index.season1"] + 1
#   status = row[-c(1)] 
#   status = status[first.index.season2:152] 
#   #print(status)
#   become.employed.index = match(2,status)
#   #print(become.employed.index)
#   status = status[become.employed.index:length(status)]
#   #print(status)
#   return.unemployed.index = match(1,status)  #index of first 0 after first 1 + index of first 1
#   #print(return.unemployed.index)
#   status = status[return.unemployed.index:length(status)]
#   #print(status)
#   eventual.employ.again.time = match(c(0,2,3),status)
#   #print(eventual.employ.again.time)
#   eventual.employ.again.time = eventual.employ.again.time[!is.na(eventual.employ.again.time)] #remove NA
#   if(length(eventual.employ.again.time)> 0){#become employed or censored before the end of the study
#     #print(min(eventual.employ.again.time))
#     return(min(eventual.employ.again.time))
#   }
#   return(NA) #check if it RETURN unemployed during season 2
# }) - 1
# df.nostatus[,"censored2"] = is.na(df.nostatus$utime2) # NA means still unemployed at the end of the study so censored.
# df.nostatus$utime2[is.na(df.nostatus$utime2)] <- 152 - (df.status$last.index.season1[is.na(df.nostatus$utime2)] + 1)# substitute NA value with the censored value (last month of follow up in season 2)
# 
# # for every row understand if the first zero was for a censoring reason (0, or 3 ) or for employment reason (2). In the first case put censored equal to true
# for(j in 1:nrow(df.status)){
#   df.nostatus[j,"censored2"] =  df.nostatus[j,"censored2"] |
#     df.status[j,df.nostatus[j,"utime2"] +1 +1] == 3 | #+1 because index inr start from 1 and our index of months starts from 0
#     df.status[j,df.nostatus[j,"utime2"] +1 +1] == 0
# }
table(df.nostatus$censored2)

#check if correct
# df.status$censored1 = df.nostatus$censored1
# df.status$censored2 = df.nostatus$censored2
# df.status$utime1 = df.nostatus$utime1
# df.status$utime2 = df.nostatus$utime2
# df.status[1,]

#create the dataset from the df.nostatus variables and add che computed utime and censored variables from df.status
data = df.nostatus

head(data)
#renames columns 
names(data) = c("id","region","education","sex","age","driving",
"trasportation","knowledge_dutch","immigrant",
                "utime1","censored1","utime2","censored2")
head(data)

#write.csv(data,"cleaned_data2samples.csv",row.names = FALSE)
library(fastDummies)
library(rockchalk)
#data = read.csv("cleaned_data2samples.csv")
#head(data)

##convert the no.status variables
data$knowledge_dutch = as.numeric(data$knowledge_dutch == 3) # boolean value if the person knows dutch (==3) or not 
data$driving = as.numeric(data$driving %in% c("B","BE")) # boolean value if the person has B driving licence
data$trasportation = as.numeric(data$trasportation== "A")# boolean value if the person has a car
data$sex = as.numeric(data$sex)-1 #boolean value if the person is male
data$censored1 = as.numeric(!data$censored1) # convert censored true in 1 and censored false in 0
data$censored2 = as.numeric(!data$censored2) 


#convert regions
region = factor(data$region)  # convert as factor
levels(region) <- list("Brabant" = c(11,12,13),
                       "Antwerp" = c(22,24),
                       "East Flanders" = c(31,34,35),
                       "West Flanders" = c(41,44,46),
                       "Limburg"= c(71,73)) # merge levels
data$region = region # save levels in the data



###education variable
data = data[!is.na(data$education),] #remove eventual rows with NA in education
data$education = as.numeric(sapply(data$education, numbers::div, m = 10) >=5) #education just as integer in base of level and than boolean if >= university degree

##education variable elaborated as suggested in pdf
# # the education variable as to be converted to boolean variables refering if the person has achive that level of education or not
# education = subset(data,select = c("id","education")) # extract id and education
# res = dummy_cols(.data = education, select_columns ="education") # create dummy variables
# level.education = apply(res, MARGIN = 1, FUN = match, x=1) #find the education level of the person
# for(i in 1:nrow(res)){
#   res[i,3:level.education[i]] = 1 # put all previous levels of education of that person as achieved (=1) (it starts from 3 becouse the first two variables are "id" and "education")
# }
# #merge the variables education with the data 
# res = res[ , -which(names(res) %in% c("id","education"))] 
# data = cbind(data, res)
# data = data[,-which(names(data) %in% c("education"))]
data$age2 = data$age + (cut.year - min.year)
head(data)
#save the cleaned data
write.csv(data,"cleaned_data_my_version2.csv",row.names = FALSE)
