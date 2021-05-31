# https://scrapingant.com/blog/scrape-dynamic-website-with-python


from bs4 import BeautifulSoup
# from scrapingant_client import ScrapingAntClient
#
# # Define URL with a dynamic web content
# url = "https://kami4ka.github.io/dynamic-website-example/"
#
# # Create a ScrapingAntClient instance
# client = ScrapingAntClient(token='c9cfbe07764442da98ecf72d6a734b55')
#
# # Get the HTML page rendered content
# page_content = client.general_request(url).content
#
# # Parse content with BeautifulSoup
# soup = BeautifulSoup(page_content)
# print(soup.find(id="test").get_text())




#
# from scrapingant_client import ScrapingAntClient
# from scrapingant_client import Cookie
#
# client = ScrapingAntClient(token='c9cfbe07764442da98ecf72d6a734b55')
#
# result = client.general_request(
#     'https://www.esrf.fr/Accelerators/Status', #''https://httpbin.org/cookies',
#     # cookies=[
#     #     Cookie(name='cookieName1', value='cookieVal1'),
#     #     Cookie(name='cookieName2', value='cookieVal2'),
#     # ]
# )
# print(result.content)
# f = open("tmp.txt","w")
# f.write(result.content)
# f.close()
# print("File tmp.txt written to disk.")
# # Response cookies is a list of Cookie objects
# # They can be used in next requests
# # response_cookies = result.cookies



import http.client

conn = http.client.HTTPSConnection("api.scrapingant.com")

headers = {
    'x-api-key': "c9cfbe07764442da98ecf72d6a734b55"
}

# conn.request("GET", "/v1/general?url=https%3A%2F%2Fexample.com", headers=headers)
conn.request("GET", "/v1/general?url=https%3A%2F%2Fwww.esrf.fr%2FAccelerators%2FStatus", headers=headers)


res = conn.getresponse()
data = res.read()

print(data.decode("utf-8"))
