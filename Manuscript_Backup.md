


---
title: "Potential Materials for Manuscript"
output: word_document
bibliography: Bib_20220818.bib
pandoc_args: ["--filter=pandoc-citeproc"]
---

## Introduction
Motivation
Estimating the tuition revenue of undergraduate students is a critical part of the budgetary planning [@trusheim2011predictive], as higher education institutions become more reliant on tuition to support operation [@barringer2016changing; @mitchell2016funding]. With tuition rates usualy being stable or stably increasing, predicting undergraduate student tuition revenue relies on acurate enrollment projection. Enrollment projection can be divided into two parts, admission of new students and retention of continuing students. In this project, we build a tool to predict weekly yield in the admission season, between end of Feburary and beginning of May. Institutional administrative units such as the Budget Office and the Admission Office can use this tool to monitor the gap between admission targets and potential outcomes and adjust the expected tuition revenue from incoming students accordingly.

Summary of previous studies
Although predicting admission yield is an important task of enrollment management in higher education institutions, there is scant literature available to describe the modeling strategy. Administrative units of the institutions often consider themselves lack of analytic capabilities and have to relie on consulting companies for the prediction task [@goenner2006predictive], but consulting companies typically view their models as proprietary and reluctant to disclose the modeling details [@desjardins2002analytic]. That being said, there were still intermissive efforts to study admission yield prediction in recent two decades [@aulck2019using; @chang2006applying; @desjardins2002analytic; @goenner2006predictive; @jamison2017applying; @sarafraz2015student;@shrestha2016offer]. Various types of predictive models are used in previous studies, including logistic regression [@chang2006applying; @desjardins2002analytic; @goenner2006predictive; @shrestha2016offer], neural network [@aulck2019using; @chang2006applying; sarafraz2015student; @shrestha2016offer] and tree-based models [@chang2006applying; @jamison2017applying]. The models were built based on training dataset and their performance was evaluated in on test dataset (out-of-sample test). The primary objective is to determine or improve the predictive accuracy of the models in the out-of-sample test, instead of investigating the effects of the input variables. 

Contribution
timing effect, covid effect, students' interest