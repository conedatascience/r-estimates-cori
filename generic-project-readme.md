# Introduction

This is the outline of an Enterprise Analytics Project. 

Feel free to use this as a template for all of your projects.

```
data analysis project
 ¦--data             
 ¦--data-raw         
 ¦--munge            
 ¦--src              
 ¦--output          
 ¦--report          
   |--report.Rmd
 ¦--README.md       
  °--project.Rproj 
```

## Automation

If you want to automate your process follow the following steps:

### Generate a DAG

A Directed Acyclic Graph or DAG will tell Airflow how you want your process to be automated.

```r
eadeploy::generate_dag(config_file = here::here("_config.yaml"))
```

### Request Automation

Your DAG will need to be put into production by the Data Science Infrastructure team. 

```r
eadeploy::request_schedule()
```
