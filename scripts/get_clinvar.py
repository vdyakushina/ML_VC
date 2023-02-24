#!/usr/bin/env python3

from ftplib import FTP
from datetime import datetime
from func import read_resources
from subprocess import run
import time, os, requests
import pandas as pd

res = os.path.dirname(os.path.realpath(__file__)) + '/resources.csv'
path = read_resources(res)
known_genes=path['known_genes']
gene_list=open(known_genes, 'r').read().strip().split('\n')

clinvar_host='ftp.ncbi.nlm.nih.gov'

def telegram_bot_sendtext(bot_message):
   bot_token = '5664464948:AAEU9xUTxP1JYHWrQ0RWbaqVM8tdyaVhCcQ'
   bot_chatID = '-1001511979986'
   send_text = 'https://api.telegram.org/bot' + bot_token + '/sendMessage?chat_id=' + bot_chatID + '&parse_mode=Markdown&text=' + bot_message
   response = requests.get(send_text)
   return response.json()

for db in ['variant_summary', 'submission_summary']:
	get_f='/pub/clinvar/tab_delimited/%s.txt.gz' % db
	db_path=path[db]
	write_d=os.path.dirname(db_path)
	write_f='%s/%s.txt.gz' % (write_d, db)
	ftp = FTP(clinvar_host)
	ftp.login()
	datef=ftp.voidcmd("MDTM %s" % get_f)
	last_changed_time = datetime.strptime(datef[4:], "%Y%m%d%H%M%S")
	ftp.close()
	current_time = datetime.fromtimestamp(os.path.getmtime(write_f))
	if last_changed_time > current_time:
		telegram_bot_sendtext("%s has updates" % db.replace('_', ' '))
		run('mv %s %s/%s_ba.txt.gz' %(write_f, write_d, db), shell=True)
		download_result=run('wget https://%s%s -O %s' % (clinvar_host,get_f, write_f), shell=True)
		if download_result.returncode==0:
			telegram_bot_sendtext("%s download success" % db.replace('_', ' '))
		else:
			telegram_bot_sendtext("%s download failed" % db.replace('_', ' '))
		if db=='variant_summary':
			t=pd.read_csv('%s' % (write_f), dtype={18:'str'}, compression='gzip', sep='\t')
			t=t[t['Assembly']=="GRCh37"]
			t=t[t['GeneSymbol'].isin(gene_list)]
			t.to_csv(db_path, compression='gzip', sep='\t', index=False)
			telegram_bot_sendtext("%s has been parsed" % db.replace('_', ' '))
	else:
		telegram_bot_sendtext("%s no updates" % db.replace('_', ' '))

