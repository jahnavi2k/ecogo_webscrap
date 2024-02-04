import azure.functions as func
import logging
import json
import os
import requests
from bs4 import BeautifulSoup
from model import *

app = func.FunctionApp(http_auth_level=func.AuthLevel.ANONYMOUS)

@app.route(route="ecogo-api/food")
def ecogo_api_food(req: func.HttpRequest) -> func.HttpResponse:
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
            'api_key': os.environ.get('api_key'),  # Replace with your API key
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
                return func.HttpResponse('API request failed: ' + str(r.status_code))
        except Exception as e:
            return func.HttpResponse('API request failed:' + str(e))
       
    else:
        return func.HttpResponse(
             "This HTTP triggered function executed successfully. Pass a name in the query string or in the request body for a personalized response.",
             status_code=200
        )



@app.route(route="ecogo-api/cosmetics")
def ecogo_api_cosmetics(req: func.HttpRequest) -> func.HttpResponse:
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



@app.route(route="ecogo-api/get-all-info")
def ecogo_api_recommend(req: func.HttpRequest) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')

    category = req.params.get("category")
    product = req.params.get("product_name")
    all_ing = req.params.get("ingredients").split("; ")
    sub_category = req.params.get("sub_category")

    if not product or not category or not all_ing or not sub_category:
        try:
            req_body = req.get_json()
        except ValueError:
            logging.error("No JSON body found when requesting information")
        else:
            category = req_body.get("category")
            product = req_body.get("product_name")
            all_ing = req_body.get("ingredients").split("; ")
            sub_category = req_body.get("sub_category")

    if product and category and all_ing and sub_category:
        try:
            if category == "cosmetics":
                sub_category = str(get_cluster(category, product))
            else:
                sub_category = str(get_cluster(category, sub_category))
            ingredients, exempts = generate_ings_exempts(all_ing)
            compounds, left = get_compounds(ingredients)
            summary = get_summary(compounds, exempts, left)
            summary['sub_category'] = sub_category
            return func.HttpResponse(json.dumps(summary))

        except Exception as e:
            return func.HttpResponse("There was an error: " + str(e))

       
    else:
        return func.HttpResponse(
             "This HTTP triggered function executed successfully. Pass the parameters in the query string or in the request body for a personalized response.",
             status_code=200
        )

@app.route(route="http_trigger", auth_level=func.AuthLevel.ANONYMOUS)
def http_trigger(req: func.HttpRequest) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')

    name = req.params.get('name')
    if not name:
        try:
            req_body = req.get_json()
        except ValueError:
            pass
        else:
            name = req_body.get('name')

    if name:
        return func.HttpResponse(f"Hello, {name}. This HTTP triggered function executed successfully.")
    else:
        return func.HttpResponse(
             "This HTTP triggered function executed successfully. Pass a name in the query string or in the request body for a personalized response.",
             status_code=200
        )