Submit-block ::= {
  contact {
    contact {
      name name {
        last "Doe",
        first "Jane",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "University of Science",
        div "Virology",
        city "Seattle",
        sub "WA",
        country "United States of America",
        street "123 University Road",
        email "jdoe@university.edu",
        postal-code "98105"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Doe",
            first "Jane",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "University of Science",
        div "Virology",
        city "Seattle",
        sub "WA",
        country "United States of America",
        street "123 University Road",
        postal-code "98105"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Doe",
              first "Jane",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "Viral Genome Sequences from University of Science"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:jdoe@university.edu"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
