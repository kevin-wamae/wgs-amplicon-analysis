### ____i. function to determine resistance profile based on count ----
# -------------------------------------------#
get_resistance_profile <- function(count, n) {
      if (is.na(count) || count == 0 && n == 0) {
        return(NA_character_)
      } else if (count == 0) {
        return("wt")
      } else if (count == 1) {
        return("single")
      } else if (count == 2) {
        return("double")
      } else if (count == 3) {
        return("triple")
      } # add more else if statements if more mutations are needed
  }


### ____ii. function to convert number to word ----
number_to_word <- function(num) {
  words <- c("zero", "single", "double", "triple", "quadruple", "quintuple",
             "sextuple", "septuple", "octuple", "nonuple", "decuple")
  return(words[num + 1])
}
