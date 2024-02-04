import azure.functions as func


def main(req: func.HttpRequest) -> func.HttpResponse:
    import logging
    import json
    import os
    import requests
    from bs4 import BeautifulSoup
    logging.info('Python HTTP trigger function processed a request.')

    product = req.params.get('product_name')
    if not product:
        try:
            req_body = req.get_json()
        except ValueError:
            logging.error("No JSON body found when requesting product name")
        else:
            product = req_body.get('product_name')

    if product:
        try:
            url = "https://incidecoder.com/"
            res = requests.get(url + "search?query={}".format(product))
            soup_data = BeautifulSoup(res.text, "html.parser")

            links = soup_data.find_all("a", class_="klavika simpletextlistitem")
            count = len(links)
          
            top_5_list={}
            
            while len(top_5_list)<5:
                i=len(top_5_list)
                if i+1 > count :
                    break
                specific = links[i]["href"]
                prod_name=links[i].text
                    # print(final_prod_name)
                    
                details = requests.get(url + specific[1:])
                soup = BeautifulSoup(details.text, "html.parser")
                if prod_name not in top_5_list:
                    top_5_list[prod_name]=[]
                
                # all_ing = ()
                for li in soup.find_all("a", class_="ingred-link black"):
                    top_5_list[prod_name].append(li.text.lower())

            return func.HttpResponse(json.dumps(top_5_list))

        except Exception as e:
            return str(e)
        
    else:
        return func.HttpResponse(
             "This HTTP triggered function executed successfully. Pass a name in the query string or in the request body for a personalized response.",
             status_code=200
        )


