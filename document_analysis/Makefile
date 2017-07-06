
all:

data:
	mkdir -p data

data/enron_mail_20150507.tgz: data
	curl --get https://www.cs.cmu.edu/~./enron/enron_mail_20150507.tgz \
	     --output data/enron_mail_20150507.tgz

data/maildir: data/enron_mail_20150507.tgz
	cd data ; tar -xkf enron_mail_20150507.tgz
