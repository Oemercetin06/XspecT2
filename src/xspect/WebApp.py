"""XspecT Flask web app"""

import re
import sys
import warnings
import subprocess
import os
import csv
import json
import time
import logging
import pickle
import secrets
import pandas as pd
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
from xspect.Classifier import classify
from xspect.search_filter import (
    single_oxa,
    get_added_genomes,
    read_search,
    read_search_spec,
    pre_processing,
    pre_processing_prefilter2,
    read_search_pre,
)
from xspect.train_filter.interface_XspecT import load_translation_dict


warnings.filterwarnings("ignore")

# Source Logging and Error Handling
# https://flask.palletsprojects.com/en/1.1.x/logging/
# https://pythonise.com/series/learning-flask/flask-error-handling
# Logging Source: https://stackoverflow.com/questions/17743019/flask-logging-cannot-get-it-to-write-to-a-file
logging.basicConfig(filename="logger.log", level=logging.ERROR)

# init WebApp with flask
app = Flask(__name__)

app.secret_key = "test"


# Initialisiere eine leere Liste für die Optionen
def load_saved_options():
    try:
        with open("saved_options.txt", "r") as file:
            options = [line.strip() for line in file]
            return options
    except FileNotFoundError:
        return []


saved_options = load_saved_options()


# Error Handling:
# https://pythonise.com/series/learning-flask/flask-error-handling


@app.route("/load_saved_options", methods=["GET"])
def load_saved_options_route():
    options = load_saved_options()
    return jsonify({"options": options})


@app.errorhandler(404)
def not_found(e):
    return render_template("404.html")


@app.errorhandler(500)
def not_found(e):
    app.logger.error(f"SERVER ERROR 500 at route {request.url} with error message: {e}")
    app.logger.error(
        f'Parameters: IC_Lookup{session.get("IC_lookup")}, \n'
        f'OXA: {session.get("OXA")}, \n'
        f'QUICK: {session.get("quick")}, \n'
        f'Filename: {session.get("filename")}, \n'
        f'Vals OXA: {session.get("vals_oxa")}, \n'
        f'Vals IC: {session.get("vals_ct")}, \n'
        f'Hits IC: {session.get("hits_ct")}, \n'
        f'Time: {session.get("time")}, \n'
        f'Prediction: {session.get("prediction")}'
    )
    return render_template("500.html")


@app.errorhandler(400)
def not_found(e):
    return render_template("400.html")


@app.errorhandler(401)
def not_found(e):
    return render_template("401.html")


# redirects to the homepage
@app.route("/", methods=["GET", "POST"])
def redirect_home():
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
    quick = session.get("quick")
    metagenome = session.get("metagenome")
    added = session.get("added", None)
    oxa = session.get("OXA", None)
    genus = session.get("genus")
    start = time.time()

    if not (os.path.exists(filename)):
        # in case that user types in route of loading screen
        # or file does not exist anymore
        return redirect("/resultsspec")

    else:
        # Checking file type
        # if the file is fasta -> concat lines
        ext = filename.split(".")[-2]
        with open(filename) as f:
            reads = f.read().splitlines()

        # Concat Lines if not .fq file
        if ext != "fq" and ext != "fastq":
            reads = "".join(reads)
            reads = reads.split(">")
            if quick:
                quick = 1
            else:
                quick = 0
            if metagenome:
                quick = 4
            reads.pop(0)
        else:
            if metagenome:
                quick = 4
            else:
                quick = 2
        # deleting file
        os.remove(filename)

        for i in range(len(reads)):
            reads[i] = reads[i].upper()
    # starts the lookup for a given sequence
    if metagenome:
        start_meta = time.time()
        reads, reads_oxa = read_search_pre(reads, BF_Master_prefilter, ext)
        end_meta = time.time()
        needed_meta = round(end_meta - start_meta, 2)
        print("Runtime filtering: ", needed_meta)

    # Lookup in Bloomfilter
    # reads should be a list of sequences

    score_ct, names_ct, hits_ct, predictions = read_search_spec(
        reads, quick, BF_Master, ext, genus
    )
    print("Testing Step 2")

    # preparing reads for oxa search
    if metagenome:
        reads = reads_oxa

    # storing values in session for creating plot
    if metagenome:
        reads_classified = score_ct
        session["reads_classified"] = reads_classified
        # assign reads to species
        species_dict = {}
        predictions_names = set()
        for ele in predictions:
            predictions_names.add(ele)
        for species in predictions_names:
            species_dict[species] = []
        # dict with species as keys and reads as values for oxa search
        for i in range(len(predictions)):
            species_dict[predictions[i]].append(reads[i])

    if not metagenome:
        session["vals_ct_spec"] = score_ct
        session["names_ct_spec"] = names_ct
        session["hits_ct_spec"] = hits_ct

    if oxa:
        if not metagenome:
            (
                score_oxa,
                names_oxa,
                coordinates_forward,
                coordinates_reversed,
            ) = single_oxa(reads, ext)
            for k in range(len(score_oxa)):
                if score_oxa[k] > 1:
                    score_oxa[k] = 1
            session["vals_oxa_spec"] = score_oxa
            session["names_oxa_spec"] = names_oxa
            session["coordinates_forward"] = coordinates_forward
            session["coordinates_reversed"] = coordinates_reversed
        elif metagenome:
            # lookup for individual species
            score_oxa_list = []
            names_oxa_list = []
            coordinates_forward_list = []
            coordinates_reversed_list = []
            for species in species_dict:
                (
                    score_oxa,
                    names_oxa,
                    coordinates_forward,
                    coordinates_reversed,
                ) = single_oxa(species_dict[species], ext)
                for k in range(len(score_oxa)):
                    if score_oxa[k] > 1:
                        score_oxa[k] = 1
                score_oxa_list.append((score_oxa, species))
                names_oxa.append(names_oxa_list)
                coordinates_forward_list.append(coordinates_forward)
                coordinates_reversed_list.append(coordinates_reversed)
            oxa_results = []
            # Ansatz mit den Tupeln geht nicht mit max funktion weil zu viele argumente
            for i in range(len(score_oxa_list)):
                if max(score_oxa_list[i][0]) > 0:
                    oxa_results.append(
                        (
                            score_oxa_list[i][1],
                            names_oxa[
                                score_oxa_list[i][0].index(max(score_oxa_list[i][0]))
                            ],
                            max(score_oxa_list[i][0]),
                        )
                    )
            # print(oxa_results)
            session["oxa_results"] = oxa_results
            session["vals_oxa_spec"] = score_oxa_list
            session["names_oxa_spec"] = names_oxa_list
            session["coordinates_forward"] = coordinates_forward_list
            session["coordinates_reversed"] = coordinates_reversed_list

    else:
        session["oxa_results"] = "None"
        session["vals_oxa_spec"] = "None"
        session["names_oxa_spec"] = "None"

    # making prediction
    if not metagenome:
        # add dynamic path
        prediction = classify(
            r"Training_data/" + genus + "_Training_data_spec.csv", score_ct, True
        )
        prediction_claast = prediction
        if prediction == "sp.":
            prediction = "NONE of the known " + genus + " species"
        else:
            translation_dict = load_translation_dict(genus)
            prediction = translation_dict[prediction]
    else:
        prediction_claast = None
        prediction = []
        for species in reads_classified:
            prediction.append(species)

    session["prediction"] = prediction

    end = time.time()
    needed = round(end - start, 2)
    print("Runtime: ", needed)
    session["time"] = str(needed)

    if prediction_claast == "470":
        IC_lookup = [True, True, True, True, True, True, True, True, False]
        score_claast, names_claast, hits_claast = read_search(IC_lookup, reads, quick=1)
        # making prediction
        prediction_claast = classify(
            r"Training_data/Training_data_IC.csv", score_claast, IC_lookup
        )
        # Making Label look nicer
        if "IC" in prediction_claast and len(prediction_claast) == 3:
            prediction_claast = "International Clone " + prediction_claast[2]
        elif prediction_claast == "None":
            prediction_claast = "NONE of the selected Clones or Genomes"
        else:
            pass
        session["prediction_claast"] = prediction_claast
        session["vals_claast"] = score_claast
        session["names_claast"] = names_claast
        session["hits_claast"] = hits_claast
        app.logger.info(
            "Assignment done for " + str(filename) + ", Time needed: " + str(needed)
        )
        return redirect("/resultsspec")
    else:
        session["prediction_claast"] = "n/a"
        session["vals_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
        session["names_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
        session["hits_claast"] = [0, 0, 0, 0, 0, 0, 0, 0]
        app.logger.info(
            "Assignment done for " + str(filename) + ", Time needed: " + str(needed)
        )
        return redirect("/resultsspec")

    app.logger.info(
        "Assignment done for " + str(filename) + ", Time needed: " + str(needed)
    )
    return redirect("/resultsspec")


# about page
@app.route("/about")
def about():
    """returns about page"""
    counter = json.load(open(r"filter/OXAs_dict/counter.txt"))
    ids = [*counter]
    r = csv.reader(open(r"Training_data/Training_data_IC.csv"))
    df = pd.DataFrame(data=list(r))
    svm_table = df.to_html(index=False, header=False)
    return render_template("about.html", svm_table=svm_table, oxa_ids=ids)


# load new BF
@app.route("/change_genus", methods=["GET", "POST"])
def change_genus():
    """Load new BF for selected genus"""
    selected_genus = request.form.get("genus")
    session["genus"] = selected_genus

    global BF_Master
    global BF_Master_prefilter
    BF_Master = pre_processing(selected_genus)
    BF_Master_prefilter = pre_processing_prefilter2(selected_genus)

    # Leere Antwort zurückgeben
    return make_response("", 200)


# train new genus
@app.route("/train_new_genus", methods=["GET", "POST"])
def train_new_genus():
    """Train new genus"""
    if request.method == "POST":
        # extract genus name from request
        genus_name = list(request.json.values())[0]

        # Sys arguments for the Python program
        system_arguments = [genus_name, "1"]

        # Run XspecT_Trainer
        subprocess.run([sys.executable, "XspecT_trainer.py"] + system_arguments)
        print("")
        print("Training done!")

        # save genus in options
        # Überprüfe, ob die Option bereits vorhanden ist
        if genus_name not in saved_options:
            print("Saving new genus: " + genus_name)
            # Füge die Option zur Liste hinzu
            saved_options.append(genus_name)

            # Speichere die Optionen in der Datei (oder Datenbank)
            with open("saved_options.txt", "a") as file:
                file.write(genus_name + "\n")

        # Erfolgreiche Antwort zurückgeben
        return redirect("/species")

    # Leere Antwort zurückgeben
    return make_response("", 200)


# species assignment page
@app.route("/species", methods=["GET", "POST"])
def species():
    """returns species page"""
    added = get_added_genomes()
    if request.method == "POST":
        data = request.json
        if data is not None:
            filename = data[-4]
            session["quick"] = data[-3]
            session["OXA"] = data[-2]
            session["metagenome"] = data[-1]
            del data[-4:]

            if not (os.path.exists("files")):
                os.mkdir("files")
            name = r"files/" + str(secrets.token_hex(8)) + filename + ".txt"

            with open(name, "w") as filehandle:
                for read in data:
                    filehandle.write("%s\n" % read)

            session["filename"] = name

            # Returning a json signal to ajax to redirect to loading page
            # the loading page then triggers the assignment process
            app.logger.info("Assignment started for " + filename)
            return json.dumps({"success": True})

        else:
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


@app.route("/resultsspec", methods=["GET", "POST"])
def resultsspec():
    """gets XspecT-Results, creates a Plot and displays them on page with further information"""

    # CALCULATING RESULTS -----------------------------------------------------

    metagenome = session.get("metagenome")

    if not metagenome:
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
        for i in range(len(values_ct)):
            dic[clonetypes[i]] = values_ct[i]
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

        # Save results in csv file
        # TODO later
    # with open(r"Results/WebApp/results_" + filename[22:-4] + ".csv", 'w', newline='') as file:
    #    file.write("XspecT Prediction, XspecT Score, ClAssT Prediction, ClAssT Score, Oxa Prediction, Oxa Score\n")
    #   for i in range(len(values_sorted)):
    #      file.write(clonetypes_sorted[i] + ", " + str(values_sorted[i]) + ", " + str(prediction_claast) + ", " + str(values_claast) + ", " + str(prediction) + ", " + str(values_oxa[i]) + "\n")

    elif metagenome:
        reads_classified = session.get("reads_classified")
        # sort reads_classified by highest value of the second element
        sorted_reads_classified = dict(
            sorted(reads_classified.items(), key=lambda x: x[1][1], reverse=True)
        )
        # get key of reads_classified with highest value of the second element from the value
        predictions = []
        values = []
        for key, value in sorted_reads_classified.items():
            predictions.append(key)
            values.append(value[1])
        clonetypes_sorted = predictions[:12]
        values_sorted = values[:12]
        prediction = predictions[0]
        maxi = values[0]
        additional_info = []
        metagenome = True
        filename = session.get("filename")

        # Save results in csv file
        # TODO later
        # with open(r"Results/WebApp/results_" + filename[22:-4] + ".csv", 'w', newline='') as file:
        #   file.write("Prediction, Score Median, Number of Contigs, Contig-Length Median, Uniqueness, Bootstrap Median\n")
        #  for key, value in sorted_reads_classified.items():
        #     file.write(key + "," + str(value[0]) + "," + str(value[1]) + "," + str(value[2]) + "," + str(value[3]) + "," + str(value[4]) + "\n")

        for key, value in sorted_reads_classified.items():
            number_of_contigs = value[1]
            value[0] = "Score Median: " + str(value[0])
            value[1] = "Number of Contigs: " + str(number_of_contigs)
            value[2] = "Contig-Length Median: " + str(value[2])
            value[3] = "Uniqueness: " + str(value[3])
            value[4] = "Bootstrap Median: " + str(value[4])
            additional_info.append(
                value[0]
                + "\n"
                + value[1]
                + "\n"
                + value[2]
                + "\n"
                + value[3]
                + "\n"
                + value[4]
            )
        text = "Detected Acinetobacter species"

        # Values of clonetypes, is None if not existing
        values_ct = session.get("vals_ct_spec")
        hits_ct = session.get("hits_ct_spec")
        clonetypes = session.get("names_ct_spec")
        values_claast = session.get("vals_claast")
        hits_claast = session.get("hits_claast")
        clonetypes_claast = session.get("names_claast")
        prediction_claast = session.get("prediction_claast")
        # Values of OXAs
        values_oxa = session.get("vals_oxa_spec")
        oxa_names = session.get("names_oxa_spec")
        oxa_results = session.get("oxa_results")

        if oxa_results != "None":
            oxa_labels = []
            oxa_data = []
            for results in oxa_results:
                oxa_labels.append("A. " + results[0] + ": " + results[1])
                oxa_data.append(results[2])
        else:
            oxa_labels = "None"
            oxa_data = "None"

    filename = session.get("filename")[22:]
    filename = os.path.splitext(filename)[0]

    # PUBMED LITERATURE SEARCH --------------------------------------------------------------------------------------------

    # Pubmed literature search Source: https://gist.github.com/bonzanini/5a4c39e4c02502a8451d
    # and https://biopython-tutorial.readthedocs.io/en/latest/notebooks/09%20-%20Accessing%20NCBIs%20Entrez%20databases.html
    Entrez.email = "xspectBIOINF@web.de"
    handle = Entrez.esearch(
        db="pubmed", sort="relevance", retmax="10", retmode="xml", term=prediction
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
        literature_content.append(paper["MedlineCitation"]["Article"]["ArticleTitle"])
        try:
            literature_abstract.append(
                paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            )
        except:
            literature_abstract.append(["No abstract available"])

    for i in range(len(literature_content)):
        literature_id.append("paper_" + str(i))

    for record in literature_info:
        literature_authors.append(record.get("AU", "?"))
        literature_journal.append(record.get("SO", "?"))

    for i in range(len(literature_authors)):
        literature_authors[i] = " ,".join(literature_authors[i])

    for i in range(len(literature_abstract)):
        literature_abstract[i] = " ".join(literature_abstract[i])

    CLEANR = re.compile("<.*?>")

    for i in range(len(literature_content)):
        literature_content[i] = re.sub(CLEANR, "", literature_content[i])
        literature_abstract[i] = re.sub(CLEANR, "", literature_abstract[i])

    literature_all = [
        literature,
        literature_content,
        literature_abstract,
        literature_authors,
        literature_journal,
        literature_id,
    ]

    if request.method == "POST":
        data = request.json
        Entrez.email = "xspectBIOINF@web.de"
        handle = Entrez.esearch(
            db="pubmed",
            sort=str(data[1]),
            retmax=str(data[0]),
            retmode="xml",
            term=prediction,
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
            literature_abstract.append(
                paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            )

        for i in range(len(literature_content)):
            literature_id.append("paper_" + str(i))

        for record in literature_info:
            literature_authors.append(record.get("AU", "?"))
            literature_journal.append(record.get("SO", "?"))

        for i in range(len(literature_authors)):
            literature_authors[i] = " ,".join(literature_authors[i])

        for i in range(len(literature_abstract)):
            literature_abstract[i] = " ".join(literature_abstract[i])

        CLEANR = re.compile("<.*?>")

        for i in range(len(literature_content)):
            literature_content[i] = re.sub(CLEANR, "", literature_content[i])
            literature_abstract[i] = re.sub(CLEANR, "", literature_abstract[i])

        literature_all = [
            literature,
            literature_content,
            literature_abstract,
            literature_authors,
            literature_journal,
            literature_id,
        ]

        return json.dumps(literature_all)

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
