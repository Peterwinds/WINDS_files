import pymysql
pymysql.install_as_MySQLdb()
from sqlalchemy import create_engine
import pandas as pd 

#=================================Connect to DB and pull data======================================#
con=create_engine('mysql://UofABEWINDS:WINDSAWSPort2020@windsdatabase-1.cdzagwevzppe.us-west-1.rds.amazonaws.com:3306/winds_test')

planting_array=pd.read_sql('SELECT * from irrigation_activity',con=con) 


planting_array.to_excel('irrigation_activity.xlsx', engine='xlsxwriter')  