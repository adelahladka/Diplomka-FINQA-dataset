#importing sas into R
#PISA dataset 2022

#CRT_SAS = Creative Thinking item data file (creative_data)
#FLT_SAS = Financial Literacy data file (flt_cog_data,flt_qqq_data, flt_tim_data)
#PISA2022 = agtegate file
#SAS_SCH_QQQ = School questionnaire data file (school_data)
#SAS_STU_QQQ = Student questionnaire data file (student_data)
#SAS_TCH_QQQ = Teacher questionnaure data file  (teacher_data)
#STU_COG_SAS = Cognite item data file (stu_cog_data)
#STU_TIM_SAS = Questionare timing data file (stu_tim_data)

install.packages('haven')
library('haven')


school_data <- read_sas("PISA2022/SAS_SCH_QQQ/SCH/cy07_msu_sch_qqq.sas7bdat", catalog_file = "PISA2022/SAS_SCH_QQQ/SCH/CY07MSU_FMT_SCH_QQQ.SAS7BCAT")
student_data <- read_sas("PISA2022/SAS_STU_QQQ/STU/cy07_msu_stu_qqq.sas7bdat", catalog_file = "PISA2022/SAS_STU_QQQ/STU/CY07MSU_FMT_STU_QQQ.SAS7BCAT")
teacher_data <- read_sas("PISA2022/SAS_TCH_QQQ/TCH/cy07_msu_tch_qqq.sas7bdat", catalog_file = "PISA2022/SAS_TCH_QQQ/TCH/CY07MSU_FMT_TCH_QQQ.SAS7BCAT")

creative_data <- read_sas("PISA2022/CRT_SAS/CY08MSP_CRT_COG.sas7bdat")

flt_cog_data <- read_sas("PISA2022/FLT_SAS/CY08MSP_FLT_COG.sas7bdat")
flt_qqq_data <- read_sas("PISA2022/FLT_SAS/CY08MSP_FLT_QQQ.SAS7BDAT")
flt_tim_data <- read_sas("PISA2022/FLT_SAS/CY08MSP_FLT_TIM.SAS7BDAT")


stu_cog_data <- read_sas("PISA2022/STU_COG_SAS/CY08MSP_STU_COG.SAS7BDAT")
stu_tim_data <- read_sas("PISA2022/STU_TIM_SAS/CY08MSP_STU_TIM.SAS7BDAT")

