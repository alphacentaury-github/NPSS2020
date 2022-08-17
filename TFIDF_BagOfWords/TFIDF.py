# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 14:11:20 2022

@author: User
"""
# import libraries 
import numpy as np
import re
import sklearn
import numpy as np
# import holoviews as hv
import nltk 
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
#from xgboost import XGBClassifier
from nltk.stem.snowball import SnowballStemmer
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
#from keras.preprocessing.text import Tokenizer
# sklean  CountVectorizer usage 
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics.pairwise import cosine_similarity
# German text stemer and stop_words 
nltk.download('punkt')
nltk.download('stopwords')
stemmer = SnowballStemmer("german")
stop_words = set(stopwords.words("german"))
#================load documents===================
#--open files 
f= open('Hitler.txt','r', encoding='UTF8')
HT_text=f.read()
f.close() 
f= open('WhiteRose.txt','r',encoding='UTF8')
WR_text=f.read()
f.close() 
# separate texts into multiple texts by seperator '#---'  
# then the first element will be a blank  
WR_docs = WR_text.split('#---')[1:]
print('number of WR texts =', len(WR_docs) )
HT_docs = HT_text.split('#---')[1:]
print('number of hitler texts =', len(HT_docs) )
all_docs = WR_docs +HT_docs # list combine 
print('number of all texts =', len(all_docs) )
#====================================================
#========How to clean and tokenize the text?==================================
def clean_text(text, for_embedding=False):
    """
        - remove any html tags (< /br> often found)
        - Keep only ASCII + European Chars and whitespace, no digits
        - remove single letter chars
        - convert all whitespaces (tabs etc.) to single wspace
        if not for embedding (but e.g. tdf-idf):
        - all lowercase
        - remove stopwords, punctuation and stemm
    """
    RE_WSPACE = re.compile(r"\s+", re.IGNORECASE)
    RE_TAGS = re.compile(r"<[^>]+>")
    RE_ASCII = re.compile(r"[^A-Za-zÀ-ž ]", re.IGNORECASE)
    RE_SINGLECHAR = re.compile(r"\b[A-Za-zÀ-ž]\b", re.IGNORECASE)
    if for_embedding:
        # Keep punctuation
        RE_ASCII = re.compile(r"[^A-Za-zÀ-ž,.!? ]", re.IGNORECASE)
        RE_SINGLECHAR = re.compile(r"\b[A-Za-zÀ-ž,.!?]\b", re.IGNORECASE)

    text = re.sub(RE_TAGS, " ", text) # remove tags
    text = re.sub(RE_ASCII, " ", text) # Keep only ASCII + European Chars and whitespace, no digits
    text = re.sub(RE_SINGLECHAR, " ", text) #remove single letter chars
    text = re.sub(RE_WSPACE, " ", text)

    word_tokens = word_tokenize(text)
    words_tokens_lower = [word.lower() for word in word_tokens]

    if for_embedding:
        # no stemming, lowering and punctuation / stop words removal
        words_filtered = word_tokens
    else:
        words_filtered = [
            stemmer.stem(word) for word in words_tokens_lower if word not in stop_words
        ]

    text_clean = " ".join(words_filtered)
    return text_clean

def print_freq_words(input_vector,vocab,n_words=10):
    """
    input vector corresponds to the frequency or tfidf values
    for each terms. 
    the pairing between word and its index is stored in vocab 
    """
    sort_index = np.argsort(-input_vector) 
    for i in range(n_words):
        idx = sort_index[i]
        print('{}:{}'.format(vocab[idx],input_vector[idx]))
#===============================================================================
#-----without prior cleaning 
vectorizer = CountVectorizer() 
corpus = all_docs 
X = vectorizer.fit_transform(corpus) # fitted vector object 
vocab = vectorizer.get_feature_names_out()
tf_mat = X.toarray() # term frequency matrix 
transformer = TfidfTransformer()
tfidf = transformer.fit_transform(tf_mat)
idf_mat = transformer.idf_
tfidf_mat = tfidf.toarray()
cosine_sim = cosine_similarity(tfidf_mat, tfidf_mat)  
print(np.round(cosine_sim,decimals=3))
print_freq_words(tf_mat[10],vocab,n_words=30) # example of Hitler's text 
  
#-----With prior cleaning 
all_docs_cleaned = [clean_text(i) for i in all_docs]
vectorizer = CountVectorizer() 
corpus = all_docs_cleaned 
X = vectorizer.fit_transform(corpus) # fitted vector object 
vocab_1 = vectorizer.get_feature_names_out()
tf_mat_1 = X.toarray() # term frequency matrix 
transformer = TfidfTransformer()
tfidf_1 = transformer.fit_transform(tf_mat_1)
idf_mat_1 = transformer.idf_
tfidf_mat_1 = tfidf_1.toarray()
cosine_sim_1 = cosine_similarity(tfidf_mat_1, tfidf_mat_1)
print(np.round(cosine_sim_1,decimals=3))
print_freq_words(tf_mat_1[10],vocab_1,n_words=30) # example of Hitler's text 

#-----if combine all docs into one file 
vectorizer = CountVectorizer()
#X = vectorizer.fit_transform([WR_text,HT_text])
X = vectorizer.fit_transform([clean_text(WR_text),clean_text(HT_text)])
tf_mat_2 = X.toarray() # term frequency matrix 
vocab_2 = vectorizer.get_feature_names_out()
print_freq_words(tf_mat_2[0],vocab_2,n_words=30) 
print_freq_words(tf_mat_2[1],vocab_2,n_words=30) 

