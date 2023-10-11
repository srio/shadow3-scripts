import urllib.request
target_url = "https://raw.githubusercontent.com/PaNOSC-ViNYL/Oasys-PaNOSC-Workspaces/master/EBS_ID32.ows"
data = urllib.request.urlopen(target_url)
for line in data:
    print(line)
