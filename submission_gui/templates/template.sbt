Submit-block ::= {
  contact {
    contact {
      name name {
        last "Zhang",
        first "Wenda",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Ningxia University",
        div "College of Life Sciences",
        city "Yinchuan",
        sub "Ningxia",
        country "China",
        street "Ningxia University",
        email "zhangwenda@example.com",
        postal-code "750021"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Zhang",
            first "Wenda",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Ningxia University",
        div "College of Life Sciences",
        city "Yinchuan",
        sub "Ningxia",
        country "China",
        street "Ningxia University",
        postal-code "750021"
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
              last "Zhang",
              first "Wenda",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "Plant virome of Lycium chinense in Ningxia, China"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:zhangwenda@example.com"
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
