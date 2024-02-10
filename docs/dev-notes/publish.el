(setq org-publish-project-alist
      '(("Technical Documentation"
         :base-directory "./"
         :publishing-function org-html-publish-to-html
         :publishing-directory "../html/"
         :section-numbers nil)))
