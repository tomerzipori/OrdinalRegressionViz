.onAttach <- function(...) {
  phrase <- random_phrase()
  packageStartupMessage(paste(strwrap(phrase), "But seriously, this small package is an experiment by me and therefore probably full of bugs", collapse = "\n"))
}


random_phrase <- function() {
  phrases <- c(
    "NP is not in P!",
    "Checkout my blog at https://tomerzipori.github.io/blog/",
    "Finger-licking!",
    "Follow the train, CJ!",
    "Totally forgot about Dre!",
    "Call you Mom!",
    "Diet Cherry Coke"
  )

  sample(phrases, 1)
}

