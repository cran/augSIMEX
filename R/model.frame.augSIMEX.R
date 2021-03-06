model.frame.augSIMEX <- function (formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
                      0L)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(stats::glm)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model
}