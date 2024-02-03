import deepchem as dc
import pubchempy as pcp
from rdkit import Chem
import requests, json
import joblib


def get_compound_url(compound_name):
    # fetch PubChem CID for the given compound name
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = json.loads(response.content)
        if "IdentifierList" in data:
            cid = data["IdentifierList"]["CID"][0]
            return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        else:
            return f"https://pubchem.ncbi.nlm.nih.gov/"
    else:
        return f"https://pubchem.ncbi.nlm.nih.gov/"

tox21_tasks = ['NR-AR',
 'NR-AR-LBD',
 'NR-AhR',
 'NR-Aromatase',
 'NR-ER',
 'NR-ER-LBD',
 'NR-PPAR-gamma',
 'SR-ARE',
 'SR-ATAD5',
 'SR-HSE',
 'SR-MMP',
 'SR-p53']


model = dc.models.GraphConvModel(len(tox21_tasks), mode='classification', model_dir='models/tox_pred')
model.restore()
cosmetic_rec = joblib.load('models/recommender.sav')
cosmetic_tfid = joblib.load('models/tfidf_vectorizer.sav')
food_rec = joblib.load('models/recommender2.sav')
food_tfid = joblib.load('models/tfidf_vectorizer2.sav')

environment_task = ['NR-AhR', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5', 'SR-HSE', 'SR-MMP', 'SR-p53']
aquatic_task = ['NR-Aromatase', 'NR-ER', 'NR-ER-LBD', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5', 'SR-HSE', 'SR-MMP']

environment_task_indices = [tox21_tasks.index(task) for task in environment_task]
aquatic_task_indices = [tox21_tasks.index(task) for task in aquatic_task]

def predict_toxicity(compound):
    featurizer = dc.feat.graph_features.ConvMolFeaturizer()
    features = featurizer([compound])
    composite_preds =  model.predict_on_batch(features)[0]
    environment_preds = composite_preds[environment_task_indices]
    aquatic_preds = composite_preds[aquatic_task_indices]
    return environment_preds, aquatic_preds, composite_preds

# def get_ingredients_from_product(product):
#    return """Water, Glycerin, Propylene glycol, Glyceryl stearate, Cetyl alcohol, Stearyl alcohol, Dimethicone, Cyclomethicone""".lower().split(", ")
# ings = get_ingredients_from_product()

def generate_ings_exempts(ings):
    list_of_ingredients = set(ings)
    exemptions = set(['water', 'aqua', 'glycerine', 'salt', 'sugar', 'malic acid', 'stearic acid', 'colour', 'dextrose', 'paprika' ,'sodium chloride', 'potassium chloride', 'magnesium chloride', 'calcium chloride', 'sodium hydroxide', 'potassium hydroxide', 'ammonium hydroxide', 'hydrochloric acid', 'sulfuric acid', 'nitric acid', 'sorbic acid',  'acetic acid', 'citric acid', 'lactic acid', 'benzoic acid', 'salicylic acid', 'urea', 'glycerin', 'propylene glycol', 'ethanol', 'isopropyl alcohol', 'hexylene glycol', 'butylene glycol', 'propanediol', 'polyethylene glycol (PEG)', 'sorbitol', 'xylitol', 'sucralose', 'saccharin', 'aspartame', 'titanium dioxide', 'iron oxide'])
    ingredients = list(list_of_ingredients - exemptions)
    exempts = list(list_of_ingredients - set(ingredients))
    return ingredients, exempts

def get_compounds(ingredients):
    compounds = {}
    left = []
    for ingredient_name in ingredients:
        try:
            result = pcp.get_compounds(ingredient_name, 'name')[0]
            compound = Chem.MolFromSmiles(result.canonical_smiles)
            compounds[compound] = ingredient_name
        except:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ingredient_name}/cids/JSON"
            response = requests.get(url)
            if response.status_code == 200:
                data = json.loads(response.content)
                if "IdentifierList" in data:
                    cid = data["IdentifierList"]["CID"][0]
                    tox_request_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/2286/JSON?cid={}".format(cid)
                    tox_response = requests.get(tox_request_url)
                    tox_data = json.loads(tox_response.text)
                    print(tox_data)
            else:
                left.append(ingredient_name)
    return compounds, left

def get_summary(compounds, exempts, left):
    summary = {}
    overall = 0
    ingredient_count = len(compounds) + len(exempts) + len(left)
    aqua_tot, env_tot = 0, 0
    for compound, name in compounds.items():
        environment_preds, aquatic_preds, composite_preds = predict_toxicity(compound)
        env, aqua, comp = environment_preds.mean(axis=0)[0], aquatic_preds.mean(axis=0)[0], composite_preds.mean(axis=0)[0]
        overall += comp
        aqua_tot += aqua
        env_tot += env
        summary[name.capitalize()] = [str(round(aqua, 3)), str(round(env, 3)), get_compound_url(name)]
    summary["Overall"] = str(round((overall / ingredient_count), 3))
    summary["Aquatic"] = str(round((aqua_tot / ingredient_count), 3))
    summary["Environment"] = str(round((env_tot / ingredient_count), 3))
    for ing in exempts:
        summary[ing.capitalize()] = ['0', '0', get_compound_url(ing)]
    for ing in left:
        summary[ing.capitalize()] = ['-1', '-1', ""]
    return summary


def get_cluster(product_type, info):
    if product_type == 'food':
        return food_rec.predict(food_tfid.transform([info]))[0]
    else:
        return cosmetic_rec.predict(cosmetic_tfid.transform([info]))[0]

