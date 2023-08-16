#!/usr/bin/env python
 
"""Script to fetch data from basespace for a given project ID.

Julien Tremblay - jtremblay514@gmail.com
"""

import os
import sys
import argparse
import re
import json
import math
import socket

#from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
from urllib2 import Request, urlopen, URLError

def restrequest(rawrequest):
   request = Request(rawrequest)

   try:
      response = urlopen(request)
      json_string = response.read()
      json_obj = json.loads(json_string)

   except URLError, e:
      print 'Got an error code:', e
      sys.exit()

   return json_obj

def downloadrestrequest(rawrequest,path):
   dirname = str(ProjectID) + os.sep + os.path.dirname(path)

   if dirname != '':
      if not os.path.isdir(dirname):
         os.makedirs(dirname)

   request = (rawrequest)

   outfile = open(str(ProjectID) + os.sep + path,'wb')
   
   try:
      response = urlopen(request,timeout=1)
      outfile.write(response.read())
      outfile.close()

   except URLError, e:
      print 'Got an error code:', e
      outfile.close()
      downloadrestrequest(rawrequest,path)


   except socket.error:
      print 'Got a socket error: retrying'
      outfile.close()
      downloadrestrequest(rawrequest,path)


def main(arguments):
 
   parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
   #parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
   #iparser.add_argument('-o', '--outfile', help="Output file", default=sys.stdout, type=argparse.FileType('w'))
   parser.add_argument('-p', '--project_id', help="Basespace Project ID", type=int)
   parser.add_argument('-r', '--run_id', help="Basespace Run ID", type=int)
   parser.add_argument('-t', '--access_token', help="Basespace access token")
    
   args = parser.parse_args(arguments)
   print args

   # initialize an authentication object using the key and secret from your app
   # Fill in with your own values
   #client_key      = '67c1c85e9cfd487c9db701d6c5e9bd6e'
   #client_secret   = 'e9279714d3c24a5b898e6ce844d75b21'
   #app_session_id  = ''
   #access_token    = '2e560b542f284938add49077bc0fb9b6'
   #basespace_url   = 'https://api.basespace.illumina.com/'
   #version         = 'v1pre3'

   #options = arg_parser()
   print args.run_id
   print args.access_token
   
   global RunID 
   RunID = str(args.run_id)
   global AccessToken 
   AccessToken = str(args.access_token)
   global ProjectID
   ProjectID = str(args.project_id)

   #request = 'http://api.basespace.illumina.com/v1pre3/runs/%s/files?access_token=%s' %(ProjectID,AccessToken)
   request = 'http://api.basespace.illumina.com/v1pre3/projects/%s/files?access_token=%s' %(ProjectID,AccessToken)
   print request

   json_obj = restrequest(request)
   totalCount = json_obj['Response']['TotalCount']

   noffsets = int(math.ceil(float(totalCount)/1000.0))

   hreflist = []
   pathlist = []
   filenamelist = [] 
   for index in range(noffsets):
      offset = 1000*index
      #request = 'http://api.basespace.illumina.com/v1pre3/runs/%s/files?access_token=%s&limit=1000&Offset=%s' %(ProjectID,AccessToken,offset) 
      request = 'http://api.basespace.illumina.com/v1pre3/projects/%s/files?access_token=%s&limit=1000&Offset=%s' %(ProjectID,AccessToken,offset) 
      json_obj = restrequest(request)

      print json_obj

      nfiles = len(json_obj['Response']['Items'])
      for fileindex in range(nfiles):
         href = json_obj['Response']['Items'][fileindex]['Href']
         hreflist.append(href)
         path = json_obj['Response']['Items'][fileindex]['Path']
         pathlist.append(path)
   
   for index in range(len(hreflist)):
      request = 'http://api.basespace.illumina.com/%s/content?access_token=%s'%(hreflist[index],AccessToken)
      print 'downloading %s' %(pathlist[index]) 
      downloadrestrequest(request, pathlist[index])
    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


