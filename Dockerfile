FROM underworldcode/underworld2:dev

CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
