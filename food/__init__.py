import azure.functions as func
import logging
import os
import json
import requests
from bs4 import BeautifulSoup

def main(req: func.HttpRequest) -> func.HttpResponse:
     
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
        # food
            url = "https://api.nal.usda.gov/fdc/v1/foods/search"
            params = {
            'api_key': os.getenv("api_key"),  # Replace with your API key
            'query': product,  # Replace with any other required parameters
            "dataType": [
                    "Branded"
                    # "Foundation",
                    # "SR Legacy"
                    ],
            "pageNumber": 1,
            "numberOfResultsPerPage": 50,
            }

            
            r = requests.get(url, params=params)
      
            if r.status_code == 200: 
                data = r.json()  
                # return data
                items = data['foods']  
                top_5_list={}
                # count=0
                foodcrit=data['foodSearchCriteria']
                count=foodcrit['numberOfResultsPerPage']
                # print(count)
                while len(top_5_list)<5 :
                    if(count<=0) :
                        break
                    count-=1
                    item=items[len(top_5_list)]
                    item_name= item['description']
                    # print(item_name)
                    if item_name not in top_5_list:
                        top_5_list[item_name]=[]
                        top_5_list[item_name].append(item['foodCategory'])
                        top_5_list[item_name].append(item['ingredients'].lower())
                
                return func.HttpResponse(json.dumps(top_5_list))             
            else:
                return func.HttpResponse('API request failed: ' + str(r.status_code), product, r)
        except Exception as e:
            return func.HttpResponse('API request failed:' + str(e))
       
    else:
        return func.HttpResponse(
             "This HTTP triggered function executed successfully. Pass a name in the query string or in the request body for a personalized response.",
             status_code=200
        )




