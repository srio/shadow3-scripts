import json
from urllib.request import urlopen

file_url = "https://raw.githubusercontent.com/srio/shadow3-scripts/master/ESRF-LIGHTSOURCES-EBS/ebs_ids.json"

u = urlopen(file_url)
ur = u.read()
url = ur.decode(encoding='UTF-8')

# print(url)

dictionary = json.loads(url)

for key in dictionary.keys():
    print(key, dictionary[key])
