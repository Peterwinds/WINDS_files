import pandas as pd
import pymysql
pymysql.install_as_MySQLdb()
from sqlalchemy import create_engine


db=create_engine('mysql://UofABEWINDS:WINDSAWSPort2020@windsdatabase-1.cdzagwevzppe.us-west-1.rds.amazonaws.com:3306/winds_test')
connection = db.raw_connection()
cursor = connection.cursor()

myfile1="Irrigation.csv"

df_file=pd.read_csv(myfile1,sep=",")


df_file.to_sql(name='irrigation_activity',con=db,if_exists='append',index=False)  #test table is the name of the table you want, it will auto append if the table exists

#if_exists{‘fail’, ‘replace’, ‘append’}, default ‘fail’
#How to behave if the table already exists.
#
#fail: Raise a ValueError.
#
#replace: Drop the table before inserting new values.
#
#append: Insert new values to the existing table.