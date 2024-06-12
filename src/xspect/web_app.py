"""XspecT Flask web app"""

# pylint: disable=function-redefined, unused-argument

from pathlib import Path
import re
import sys
import warnings
import subprocess
import os
import json
import time
import logging
import secrets
from Bio import Entrez, Medline
from flask import (
    Flask,
    render_template,
    session,
    request,
    redirect,
    abort,
    make_response,
    jsonify,
)
import xspect.model_management as mm
from xspect.definitions import get_xspect_upload_path
from xspect.train import train_ncbi


warnings.filterwarnings("ignore")

# Source Logging and Error Handling
logging.basicConfig(filename="logger.log", level=logging.ERROR)

# init WebApp with flask
app = Flask(__name__)

app.secret_key = "test"


saved_options = mm.get_models()["Species"] if mm.get_models() else []


# Error Handling:
# https://pythonise.com/series/learning-flask/flask-error-handling


@app.route("/load_saved_options", methods=["GET"])
def load_saved_options_route():
    """Loads saved options"""
    options = mm.get_models()["Species"]
    return jsonify({"options": options})


@app.errorhandler(404)
def not_found(e):
    """Error handling for 404"""
    return render_template("404.html")


@app.errorhandler(500)
def not_found(e):
    """Error handling for 500"""
    app.logger.error(
        "SERVER ERROR 500 at route {url} with error message: {error}",
        url=request.url,
        error=e,
    )
    return render_template("500.html")


@app.errorhandler(400)
def not_found(e):
    """Error handling for 400"""
    return render_template("400.html")


@app.errorhandler(401)
def not_found(e):
    """Error handling for 401"""
    return render_template("401.html")


@app.route("/", methods=["GET", "POST"])
def redirect_home():
    """Redirects to the homepage"""
    return redirect("home")


# about page
@app.route("/home")
def home():
    """returns home page"""
    return render_template("home.html")


# Starts Assignment-Process for AspecT and leads to result-page
@app.route("/assignspec")
def assignspec():
    """Uses User Options to process the file, returns a signal to the loadingpage to go the the
    result-page when done"""

    # getting user parameters back with session function
    filename = session.get("filename", None)
    metagenome = session.get("metagenome")
    genus = session.get("genus")
    start = time.time()

    if not os.path.exists(filename):
        # in case that user types in route of loading screen
        # or file does not exist anymore
        return redirect("/resultsspec")

    sequence_input = Path(filename)
    species_filter_model = mm.get_species_model(genus)

    if metagenome:
        genus_filter_model = mm.get_genus_model(genus)
        filtered_sequences = genus_filter_model.filter(sequence_input)
        prediction, scores = species_filter_model.predict(
            filtered_sequences["Acinetobacter"]
        )
    else:
        prediction, scores = species_filter_model.predict(sequence_input)

    prediction = species_filter_model.display_names[prediction[0]]

    session["vals_ct_spec"] = list(scores.values())
    session["names_ct_spec"] = [
        species_filter_model.display_names[i] for i in scores.keys()
    ]
    session["hits_ct_spec"] = "None"

    session["oxa_results"] = "None"
    session["vals_oxa_spec"] = "None"
    session["names_oxa_spec"] = "None"

    session["prediction"] = prediction

    end = time.time()
    needed = round(end - start, 2)
    print("Runtime: ", needed)
    session["time"] = str(needed)

    session["prediction_claast"] = "n/a"
    session["vals_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
    session["names_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
    session["hits_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
    app.logger.info(
        "Assignment done for {file}, Time needed: {time}",
        file=str(filename),
        time=str(needed),
    )
    return redirect("/resultsspec")


# about page
@app.route("/about")
def about():
    """returns about page"""
    return render_template("about.html", svm_table=None, oxa_ids=None)


# load new BF
@app.route("/change_genus", methods=["GET", "POST"])
def change_genus():
    """Load new BF for selected genus"""
    selected_genus = request.form.get("genus")
    session["genus"] = selected_genus
    return make_response("", 200)


# train new genus
@app.route("/train_new_genus", methods=["GET", "POST"])
def train_new_genus():
    """Train new genus"""
    if request.method == "POST":
        # extract genus name from request
        genus_name = list(request.json.values())[0]

        # Run XspecT_Trainer
        train_ncbi(genus_name, "1", "", "")
        print("")
        print("Training done!")

        # save genus in options
        # Überprüfe, ob die Option bereits vorhanden ist
        if genus_name not in saved_options:
            print("Saving new genus: " + genus_name)
            # Füge die Option zur Liste hinzu
            saved_options.append(genus_name)

        # Erfolgreiche Antwort zurückgeben
        return redirect("/species")

    # Leere Antwort zurückgeben
    return make_response("", 200)


# species assignment page
@app.route("/species", methods=["GET", "POST"])
def species():
    """returns species page"""
    added = None
    if request.method == "POST":
        data = request.json
        if data is not None:
            filename = data[-4]
            session["quick"] = data[-3]
            session["OXA"] = data[-2]
            session["metagenome"] = data[-1]
            del data[-4:]

            upload_path = get_xspect_upload_path() / (
                str(secrets.token_hex(8)) + filename
            )

            with open(upload_path, "w", encoding="utf-8") as filehandle:
                for read in data:
                    filehandle.write(f"{read}\n")

            session["filename"] = str(upload_path)

            # Returning a json signal to ajax to redirect to loading page
            # the loading page then triggers the assignment process
            app.logger.info("Assignment started for {file}", file=filename)
            return json.dumps({"success": True})

        # Source: https://flask-restplus.readthedocs.io/en/stable/errors.html
        abort(400)
    return render_template(
        "species.html",
        added=added,
        results_oxa=[0, 0, 0, 0],
        oxas="None",
        results_ct=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        hits_ct=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        clonetypes=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        results_claast=[0, 0, 0, 0, 0, 0, 0, 0],
        hits_claast=[0, 0, 0, 0, 0, 0, 0, 0],
        clonetypes_claast=[0, 0, 0, 0, 0, 0, 0, 0],
        filename="filename",
        maxi=1,
        time=0,
        prediction="n/a",
        prediction_claast="n/a",
        literature="",
        literature_content="",
        literature_abstract="",
        literature_authors=[[""], [""], [""], [""], [""], [""], [""], [""], [""], [""]],
        literature_journal="",
        literature_all="",
        text="",
        additional_info="",
        metagenome=False,
        oxa_labels="",
        oxa_data="",
    )


def literature_search(query):
    """Searches for literature on PubMed"""
    # Pubmed literature search Source: https://gist.github.com/bonzanini/5a4c39e4c02502a8451d
    Entrez.email = "xspectBIOINF@web.de"
    try:
        handle = Entrez.esearch(
            db="pubmed", sort="relevance", retmax="10", retmode="xml", term=query
        )
        pubmed_results = Entrez.read(handle)

        id_list = pubmed_results["IdList"]
        literature = []
        for i in id_list:
            literature.append("https://pubmed.ncbi.nlm.nih.gov/" + str(i) + "/")
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
        papers = Entrez.read(handle)

        handle2 = Entrez.efetch(db="pubmed", id=ids, rettype="medline")
        literature_info = Medline.parse(handle2)
        literature_info = list(literature_info)

        literature_content = []
        literature_abstract = []
        literature_authors = []
        literature_journal = []
        literature_id = []
        for paper in papers["PubmedArticle"]:
            literature_content.append(
                paper["MedlineCitation"]["Article"]["ArticleTitle"]
            )
            try:
                literature_abstract.append(
                    paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
                )
            except KeyError:
                literature_abstract.append(["No abstract available"])

        for i in range(len(literature_content)):
            literature_id.append("paper_" + str(i))

        for record in literature_info:
            literature_authors.append(record.get("AU", "?"))
            literature_journal.append(record.get("SO", "?"))

        literature_authors = [" ,".join(authors) for authors in literature_authors]
        literature_abstract = [" ".join(abstract) for abstract in literature_abstract]

        cleanr = re.compile("<.*?>")
        literature_content = [re.sub(cleanr, "", i) for i in literature_content]
        literature_abstract = [re.sub(cleanr, "", i) for i in literature_abstract]

        literature_all = [
            literature,
            literature_content,
            literature_abstract,
            literature_authors,
            literature_journal,
            literature_id,
        ]

        return literature_all
    except Exception:
        return None


@app.route("/resultsspec", methods=["GET", "POST"])
def resultsspec():
    """gets XspecT-Results, creates a Plot and displays them on page with further information"""

    # Values of clonetypes, is None if not existing
    filename = session.get("filename")
    values_ct = session.get("vals_ct_spec")
    hits_ct = session.get("hits_ct_spec")
    clonetypes = session.get("names_ct_spec")
    values_claast = session.get("vals_claast")
    hits_claast = session.get("hits_claast")
    clonetypes_claast = session.get("names_claast")
    prediction = session.get("prediction")
    prediction_claast = session.get("prediction_claast")
    # Values of OXAs
    values_oxa = session.get("vals_oxa_spec")
    oxa_names = session.get("names_oxa_spec")
    additional_info = "Score"
    maxi = 1
    text = "Most similar Acinetobacter species"
    metagenome = False
    oxa_labels = "None"
    oxa_data = "None"

    dic = {}
    clonetypes_sorted = []
    # the values will be sorted by highest values for better readability
    for idx, ct_value in enumerate(values_ct):
        dic[clonetypes[idx]] = ct_value
    values_sorted = sorted(values_ct, reverse=True)
    for i in sorted(dic, key=dic.get, reverse=True):
        clonetypes_sorted.append(i)

    # only the 10 biggest values will be shown for visibility
    if len(values_sorted) > 10:
        values_sorted = values_sorted[:10]
        clonetypes_sorted = clonetypes_sorted[:10]

    # if less then 5 values are found, add empty values
    if len(values_sorted) < 5:
        for i in range(5 - len(values_sorted)):
            values_sorted.append(0)
            clonetypes_sorted.append("n/a")

    filename = session.get("filename")[22:]
    filename = os.path.splitext(filename)[0]

    literature_all = literature_search(prediction)

    return render_template(
        "species.html",
        results_oxa=values_oxa,
        oxas=oxa_names,
        results_ct=values_sorted,
        hits_ct=hits_ct,
        clonetypes=clonetypes_sorted,
        results_claast=values_claast,
        hits_claast=hits_claast,
        clonetypes_claast=clonetypes_claast,
        filename=filename,
        maxi=maxi,
        time=session.get("time"),
        prediction=prediction,
        prediction_claast=prediction_claast,
        literature_all=literature_all,
        additional_info=additional_info,
        text=text,
        metagenome=metagenome,
        oxa_labels=oxa_labels,
        oxa_data=oxa_data,
    )
