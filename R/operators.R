#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Tee pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%T>\%}} for details.
#'
#' @name %T>%
#' @rdname tee
#' @keywords internal
#' @export
#' @importFrom magrittr %T>%
#' @usage lhs \%T>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The left-hand side, `lhs`.
#' @details The tee pipe works like %>%, except the return value is lhs itself, and not the result of rhs function/expression.
NULL

#' Exposition pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%$\%}} for details.
#'
#' @name %$%
#' @rdname exposition
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
#' @param lhs A list, environment, or a data.frame.
#' @param rhs An expression where the names in lhs is available.
#' @return The result of evaluating `rhs` with the contents of `lhs` exposed.
#' @details This operator exposes the contents of the left-hand side object to the expression on the right side.
NULL

#' Not %in% operator
#'
#' See \code{purrr::\link[purrr:negate]{negate}} for details.
#'
#' @name %ni%
#' @rdname ni-operator
#' @keywords internal
#' @export
#' @importFrom purrr negate
#' @return A logical vector, `TRUE` where an element of the left-hand side is not in the right-hand side.
#' @details The opposite of %in%
'%ni%' <- purrr::negate(`%in%`)
