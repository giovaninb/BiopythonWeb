#!/usr/bin/env python
# -*- coding: utf-8 -*-

# created by giggio in 05/12/2017

import os
from os.path import join, dirname, realpath
from flask import (url_for, flash, redirect, render_template,
                   request, send_from_directory)
from werkzeug import secure_filename
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from app import app

UPLOAD_FOLDER = '/home/giovani/BiopythonWeb/app/uploads'

ALLOWED_EXTENSIONS = set(
    ['txt', 'pdf', 'fasta']
)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route("/index")
@app.route("/")
def index():
    return render_template('index.html')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':

        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)

        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return handle_fa(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    return render_template('index.html')

def handle_fa(filename):
    fa = SeqIO.parse(filename, 'fasta')
    for i in fa:
        dna = str(i.seq)
        rna = Seq(dna, IUPAC.unambiguous_dna)
        proteina = rna.translate()
        print('Sequência de DNA : %s' % rna )
        print('Proteina: %s' % proteina)
    return render_template('proteina.html', proteina=proteina, dnaSequence=rna)
