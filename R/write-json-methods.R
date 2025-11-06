# Copyright (C) 2025, Alvaro Uzaheta - SNlab-ETH Zurich
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT
# License for more details.
#
# You should have received a copy of the MIT License along with this
# program. If not, see <https://opensource.org/licenses/MIT>.

# JSON Export Helper Functions -----------------------------------------------

#' Transform a object to a valid JSON text
#'
#' @param object the object to transform.
#' @param first_position the position in the string until where the text
#'   is replaced.
#' @param first_replace the string to replace the text until the
#'   `first_position`.
#' @param last_position the position from where the text is replaced
#'   until the end of the string.
#' @param last_replace the string to replace the text from the
#'   `last_position` until the end of the string.
#' @return a valid JSON text.
#' @noRd
transform_json <- function(
  object,
  first_position = 1,
  first_replace = "",
  last_position = -1,
  last_replace = ",",
  pretty = TRUE,
  digits = NA,
  auto_unbox = TRUE,
  always_decimal = FALSE,
  factor = "integer"
) {
  json_text <- jsonlite::toJSON(
    object,
    pretty = pretty,
    digits = digits,
    auto_unbox = auto_unbox,
    always_decimal = always_decimal,
    factor = factor
  )
  stringr::str_sub(json_text, 1, first_position) <- first_replace
  stringr::str_sub(json_text, last_position) <- last_replace
  json_text
}

create_json_chunk_matrix <- function(
  indices,
  matrix,
  first_position = 3,
  first_replace = "",
  last_position = -2,
  last_replace = ","
) {
  transform_json(
    matrix[indices, , drop = FALSE],
    first_position = first_position,
    first_replace = first_replace,
    last_position = last_position,
    last_replace = last_replace
  )
}

create_json_chunk_vector <- function(
  indices,
  vector,
  first_position = 1,
  first_replace = "",
  last_position = -1,
  last_replace = ","
) {
  transform_json(
    vector[indices],
    first_position = first_position,
    first_replace = first_replace,
    last_position = last_position,
    last_replace = last_replace
  )
}

approx_nchar_vector <- function(n, last_val, indent = 0L, is_integer = TRUE) {
  avg_digits <- max(ceiling(log10(last_val)), 1)
  extra <- if (is_integer) 0 else 6
  per_element_chars <- indent + avg_digits + extra + 2 # "  123,\n"

  n * per_element_chars
}

approx_nchar_matrix <- function(
  n,
  p,
  max_val,
  is_integer = TRUE,
  indent = 0L
) {
  digit_width <- nchar(abs(max_val))
  extra <- if (is_integer) 0 else 6
  number_width <- digit_width + extra
  per_number <- indent + number_width + 2
  row_overhead <- indent + 4 # [\n, indent, ], \n
  outer_overhead <- 4
  n * p * per_number + n * row_overhead + outer_overhead
}

# S3 Generic Function ---------------------------------------------------------

#' Export a Stan data slot to JSON format
#'
#' Generic function to handle different data types when exporting Stan data
#' to JSON format. Automatically chooses the appropriate serialization method
#' based on data type and size.
#'
#' @param data The data object to export
#' @param slot_name Character string with the name of the data slot
#' @param n_chunks Number of chunks for large data (passed from write_json)
#' @param size_threshold Size threshold for chunking (2^31 - 1 by default)
#' @param ... Additional arguments passed to methods
#'
#' @return A character vector with JSON text representation
#' @export
export_stan_data_slot <- function(
  data,
  slot_name,
  n_chunks = 10,
  size_threshold = 2^31 - 1,
  ...
) {
  UseMethod("export_stan_data_slot")
}

# S3 Methods ------------------------------------------------------------------

#' @export
export_stan_data_slot.numeric <- function(
  data,
  slot_name,
  n_chunks = 10,
  size_threshold = 2^31 - 1,
  ...
) {
  is_integer <- rlang::is_integerish(data)

  data_length <- length(data)

  cli::cli_progress_step(
    c(
      "Exporting vector {.field {slot_name}}",
      "({prettyunits::pretty_num(data_length)})"
    )
  )

  approx_size <- approx_nchar_vector(
    data_length,
    max(abs(data), na.rm = TRUE),
    is_integer = is_integer
  )

  if (approx_size > size_threshold) {
    while (approx_size / n_chunks > size_threshold) {
      n_chunks <- n_chunks * 2
    }
    cli::cli_progress_message(
      "Vector {.field {slot_name}} is large, using {n_chunks} chunks"
    )
    vector_partition <- parallel::splitIndices(data_length, n_chunks)
    vector_text <- c(
      glue::glue('"{slot_name}": ['),
      vapply(
        vector_partition,
        create_json_chunk_vector,
        character(1),
        vector = data
      ),
      "],"
    )
    stringr::str_sub(vector_text[n_chunks + 1], -1) <- ""
    return(vector_text)
  } else {
    return(c(
      glue::glue('"{slot_name}": ['),
      transform_json(
        data,
        first_position = 1,
        first_replace = "",
        last_position = -1,
        last_replace = "],"
      )
    ))
  }
}

#' @export
export_stan_data_slot.list <- function(
  data,
  ...
) {
  names_slots <- cli::cli_vec(names(data), list("vec-trunc" = 3))
  cli::cli_progress_step(
    "Exporting list with slots: {.field {names_slots}}"
  )
  transform_json(
    data,
    first_position = 3,
    first_replace = "",
    last_position = -2,
    last_replace = ""
  )
}

#' @export
export_stan_data_slot.matrix <- function(
  data,
  slot_name,
  n_chunks = 10,
  size_threshold = 2^31 - 1,
  ...
) {
  dim_data <- dim(data)
  n_rows <- dim_data[1]
  n_cols <- dim_data[2]
  is_integer <- rlang::is_integerish(data)

  cli::cli_progress_step(
    c(
      "Exporting matrix {.field {slot_name}}",
      "({prettyunits::pretty_num(n_rows)} x {prettyunits::pretty_num(n_cols)})"
    )
  )

  approx_size <- approx_nchar_matrix(
    n_rows,
    n_cols,
    max_val = max(abs(data)),
    is_integer = is_integer
  )

  if (approx_size > size_threshold) {
    while (approx_size / n_chunks > size_threshold) {
      n_chunks <- n_chunks * 2
    }
    cli::cli_progress_message(
      "Matrix {.field {slot_name}} is large, using {n_chunks} chunks"
    )
    row_partition <- parallel::splitIndices(n_rows, n_chunks)
    matrix_text <- c(
      glue::glue('"{slot_name}": ['),
      vapply(
        row_partition,
        create_json_chunk_matrix,
        character(1),
        matrix = data
      ),
      "],"
    )
    stringr::str_sub(matrix_text[n_chunks + 1], -1) <- ""
    return(matrix_text)
  } else {
    return(c(
      glue::glue('"{slot_name}": ['),
      transform_json(
        data,
        first_position = 2,
        first_replace = "",
        last_position = -1,
        last_replace = "],"
      )
    ))
  }
}


#' @export
export_stan_data_slot.default <- function(
  data,
  slot_name,
  n_chunks = 10,
  size_threshold = 2^31 - 1,
  ...
) {
  if (is.null(data)) {
    return("")
  }
  cli::cli_abort(
    "There is not an exporting method for class {.cls {class(data)}} {.field {slot_name}}"
  )
}

# Enhanced write_json function -----------------------------------------------

#' Write model data to a JSON file
#'
#' When data objects are too large, their JSON representation may exceed
#' R's character string limit (2^31-1 bytes), causing write failures.
#' This function prevents such errors by dividing data into smaller chunks.
#'
#' @details
#' The chunking process operates as follows:
#' \enumerate{
#'   \item A size is estimated for each slot in the `data_stan` object
#'   \item If a slot's estimated JSON size exceeds `size_threshold`, it is
#'         divided into chunks
#'   \item The number of chunks follows the sequence: n_chunks, 2*n_chunks,
#'         4*n_chunks, 8*n_chunks, etc. (powers of 2)
#'   \item Division continues until each chunk's estimated size falls below
#'         the threshold
#'    \item The data or their chunks are converted to JSON format
#' }
#'
#' @param data_stan An object containing data slots to be serialized.
#' @param file_name File name to write the JSON string to,
#'   needs the .json extension.
#' @param n_chunks Initial number of chunks to attempt. Default is 10.
#' @param grain_size Number of events to process in each thread when
#'   using within chain parallelization with Map-Reduce.
#' @param exclude_slots Character vector of slot names to exclude from export.
#'   Default to empty character vector.
#' @param size_threshold Maximum estimated size (in characters) allowed
#'   per chunk. Default is typically set to avoid the 2^31-1 character limit.
#' @return a charater string with the file name of the JSON file
#' @export
write_json <- function(
  data_stan,
  file_name,
  n_chunks = 10,
  grain_size = 10,
  exclude_slots = character(0),
  size_threshold = 2^31 - 1
) {
  data_stan <- data_stan$data_stan
  # Add grain_size to data
  data_stan[["grain_size"]] <- grain_size

  # Get all slot names, excluding specified ones
  all_slots <- setdiff(names(data_stan), exclude_slots)

  # Split between length one slots and others
  lenght_slots <- sapply(data_stan, length)

  one_slots <- intersect(names(data_stan)[lenght_slots == 1], all_slots)
  more_slots <- setdiff(all_slots, one_slots)

  # ToDo: modify to use as default cmdstanr::write_stan_json()
  #   when the estimated size is smaller that 2^31 -1

  cli::cli_progress_bar(
    "Exporting Stan data slots",
    total = length(more_slots) + 1
  )

  # Write all parts to file
  conn <- file(file_name, open = "wb")
  on.exit(close(conn))
  writeLines("{", conn, useBytes = TRUE)

  # Export each slot using appropriate S3 method
  for (slot_name in more_slots) {
    cli::cli_progress_update()
    export_stan_data_slot(
      data = data_stan[[slot_name]],
      slot_name = slot_name,
      n_chunks = n_chunks,
      size_threshold = size_threshold
    ) |>
      writeLines(conn, useBytes = TRUE)
  }

  # Export together slots of length one
  cli::cli_progress_update()
  export_stan_data_slot.list(
    data = data_stan[one_slots]
  ) |>
    writeLines(conn, useBytes = TRUE)

  cli::cli_progress_done()

  writeLines("}", conn)

  cli::cli_inform(
    "Successfully exported {length(all_slots)} slots to {.file {file_name}}"
  )

  file_name
}

#' Write model data to a JSON file, legacy function
#'
#' When the data is too large, the JSON string can be too long and cause
#' and error when writing the data to a file
#' (R character strings are limited to 2^31-1 bytes).
#' In that case, the JSON string is divided into chunks.
#' @param x Data object
#' @param file_name File name to write the JSON string to,
#'   needs the .json extension.
#' @param n_chunks Number of chunks to divide the JSON string to avoid
#'   the error.
#' @param grain_size Number of events to process in each thread when
#'   using within chain parallelization with Map-Reduce.
#' @param export_all_slots Logical, whether to export all slots in data_stan
#'   (TRUE) or use legacy behavior (FALSE, default for backward compatibility).
#' @param exclude_slots Character vector of slot names to exclude from export
#'   when export_all_slots = TRUE.
#' @return NULL
#' @noRd
write_json_legacy <- function(
  x,
  file_name,
  n_chunks = 10,
  grain_size = 5,
  export_all_slots = FALSE,
  exclude_slots = character(0)
) {
  model <- attr(x, "model")
  sub_model <- attr(x, "sub_model")
  has_sample <- attr(x, "sample")
  scale <- attr(x, "scale")

  # object with information for postprocessing
  data_gathered <- x[!grepl("data_stan", names(x))]
  data_gathered$json_file <- file_name

  # Preserve class and attributes from the original data object
  class(data_gathered) <- class(x)
  attr(data_gathered, "model") <- model
  attr(data_gathered, "subModel") <- sub_model
  attr(data_gathered, "sample") <- has_sample
  attr(data_gathered, "scale") <- scale
  attr(data_gathered, "json_file") <- TRUE

  # prepare data for JSON and calculate approx. size
  data_stan <- x$data_stan
  n_size <- data_stan[[glue("N_{sub_model}")]]
  t_size <- data_stan[[glue("T_{sub_model}")]]
  p_size <- data_stan[[glue("P_{sub_model}")]]
  q_size <- data_stan[[glue("Q_{sub_model}")]]
  x_name <- glue("X_{sub_model}")
  z_name <- glue("Z_{sub_model}")
  has_z <- q_size > 0
  has_interaction <- !is.null(data_stan[["interaction"]])
  is_integer_x <- rlang::is_integerish(data_stan[[x_name]])

  approx_size_x <- approx_nchar_matrix(
    n_size,
    p_size,
    max_val = max(data_stan[[x_name]]),
    is_integer = is_integer_x
  )

  if (has_z) {
    is_integer_z <- rlang::is_integerish(data_stan[[z_name]])
    approx_size_z <- approx_nchar_matrix(
      n_size,
      q_size,
      max_val = max(data_stan[[z_name]]),
      is_integer = is_integer_z
    )
  } else {
    approx_size_z <- 0
    data_stan[[z_name]] <- NULL
  }
  if (has_interaction) {
    approx_size_interaction <- approx_nchar_vector(t_size, 100)
  } else {
    approx_size_interaction <- 0
  }
  if (sub_model == "rate") {
    mode(data_stan[["is_dependent"]]) <- "integer"
    approx_size_is_dependent <- approx_nchar_vector(t_size, 1)
    is_integer_timespan <- rlang::is_integerish(data_stan[["timespan"]])
    approx_size_timespan <- approx_nchar_vector(
      t_size,
      max(data_stan[["timespan"]]),
      is_integer = is_integer_timespan
    )
    approx_size_rate <- approx_size_is_dependent + approx_size_timespan
  } else {
    approx_size_rate <- 0
  }
  approx_max_size <- max(approx_size_x, approx_size_z)
  aprox_size <- approx_nchar_vector(t_size, n_size)
  total_size <- aprox_size *
    3 +
    approx_size_x +
    approx_size_z +
    approx_size_rate +
    approx_size_interaction +
    200
  if (total_size < (2^31 - 1)) {
    data_stan[["grain_size"]] <- grain_size
    cmdstanr::write_stan_json(
      data = data_stan,
      file = file_name
    )
    return(data_gathered)
  }

  if (approx_max_size > (2^31 - 1)) {
    row_partition <- parallel::splitIndices(n_size, n_chunks)
    x_text <- c(
      glue("\"X_{sub_model}\": ["),
      vapply(
        row_partition,
        create_json_chunk_matrix,
        character(1),
        matrix = data_stan[[x_name]]
      ),
      "],"
    )
    stringr::str_sub(x_text[n_chunks + 1], -1) <- ""
    if (has_z) {
      z_text <- c(
        glue("\"Z_{sub_model}\": ["),
        vapply(
          row_partition,
          create_json_chunk_matrix,
          character(1),
          matrix = data_stan[[z_name]]
        ),
        "],"
      )
      stringr::str_sub(z_text[n_chunks + 1], -1) <- ""
    }
  } else {
    x_text <- c(
      glue("\"X_{sub_model}\":"),
      toJSON(data_stan[[x_name]], pretty = TRUE, digits = NA),
      ","
    )
    if (has_z) {
      z_text <- c(
        glue("\"Z_{sub_model}\":"),
        toJSON(data_stan[[z_name]], pretty = TRUE, digits = NA),
        ","
      )
    }
  }

  if (aprox_size > (2^31 - 1)) {
    vector_partition <- parallel::splitIndices(t_size, n_chunks)
    start_text <- c(
      glue("\"start_{sub_model}\": ["),
      vapply(
        vector_partition,
        create_json_chunk_vector,
        character(1),
        vector = data_stan[[glue("start_{sub_model}")]]
      ),
      "],"
    )
    stringr::str_sub(start_text[n_chunks + 1], -1) <- ""
    end_text <- c(
      glue("\"end_{sub_model}\": ["),
      vapply(
        vector_partition,
        create_json_chunk_vector,
        character(1),
        vector = data_stan[[glue("end_{sub_model}")]]
      ),
      "],"
    )
    stringr::str_sub(end_text[n_chunks + 1], -1) <- ""
    chose_text <- c(
      glue("\"chose_{sub_model}\": ["),
      vapply(
        vector_partition,
        create_json_chunk_vector,
        character(1),
        vector = data_stan[[glue("chose_{sub_model}")]]
      ),
      "],"
    )
    stringr::str_sub(chose_text[n_chunks + 1], -1) <- ""
    sender_text <- c(
      "\"sender\": [",
      vapply(
        row_partition,
        create_json_chunk_vector,
        character(1),
        vector = data_stan[["sender"]]
      ),
      "],"
    )
    stringr::str_sub(sender_text[n_chunks + 1], -1) <- ""
    # start_group_text <- c(
    #   "\"start_group\": [",
    #   vapply(
    #     row_partition,
    #     create_json_chunk_vector,
    #     character(1),
    #   vector = data_stan[["start_group"]]
    # ),
    # "],"
    # )
    # stringr::str_sub(start_group_text[n_chunks + 1], -1) <- ""
    if (has_interaction) {
      interaction_text <- c(
        "\"interaction\": [",
        vapply(
          vector_partition,
          create_json_chunk_vector,
          character(1),
          vector = data_stan[["interaction"]]
        ),
        ","
      )
      stringr::str_sub(interaction_text[n_chunks + 1], -1) <- ""
    }
    if (sub_model == "rate") {
      timespan_text <- c(
        "\"timespan\": [",
        vapply(
          vector_partition,
          create_json_chunk_vector,
          character(1),
          vector = data_stan[["timespan"]]
        ),
        "],"
      )
      stringr::str_sub(timespan_text[n_chunks + 1], -1) <- ""
      is_dependent_text <- c(
        "\"is_dependent\": [",
        vapply(
          vector_partition,
          create_json_chunk_vector,
          character(1),
          vector = data_stan[["is_dependent"]]
        ),
        "],"
      )
      stringr::str_sub(is_dependent_text[n_chunks + 1], -1) <- ""
    }
  } else {
    start_text <- c(
      glue("\"start_{sub_model}\":"),
      toJSON(data_stan[[glue("start_{sub_model}")]], pretty = TRUE),
      ","
    )
    end_text <- c(
      glue("\"end_{sub_model}\":"),
      toJSON(data_stan[[glue("end_{sub_model}")]], pretty = TRUE),
      ","
    )
    chose_text <- c(
      glue("\"chose_{sub_model}\":"),
      toJSON(data_stan[[glue("chose_{sub_model}")]], pretty = TRUE),
      ","
    )
    sender_text <- c(
      "\"sender\":",
      toJSON(data_stan[["sender"]], pretty = TRUE),
      ","
    )
    # start_group_text <- c(
    #   "\"start_group\":",
    #   toJSON(data_stan[["start_group"]], pretty = TRUE),
    #   ","
    # )
    if (has_interaction) {
      interaction_text <- c(
        "\"interaction\":",
        toJSON(data_stan[["interaction"]], pretty = TRUE),
        ","
      )
    }
    if (sub_model == "rate") {
      timespan_text <- c(
        "\"timespan\":",
        toJSON(data_stan[["timespan"]], pretty = TRUE),
        ","
      )
      is_dependent_text <- c(
        "\"is_dependent\":",
        toJSON(data_stan[["is_dependent"]], pretty = TRUE),
        ","
      )
    }
  }

  keep_dttxt <- c(
    glue("{data}_{sub_model}", data = c("N", "T", "P", "Q")),
    c("A", if (has_interaction) "C" else NULL)
  )
  data_text <- transform_json(
    c(data_stan[keep_dttxt], list(grain_size = grain_size)),
    first_position = 3,
    first_replace = "",
    last_replace = "}"
  )

  # Write the data to a json file
  conn <- file(file_name, open = "wb")
  writeLines("{", conn, useBytes = TRUE)
  writeLines(start_text, conn)
  writeLines(end_text, conn)
  writeLines(chose_text, conn)
  writeLines(sender_text, conn)
  # writeLines(start_group_text, conn)
  writeLines(x_text, conn)
  if (has_z) {
    writeLines(z_text, conn)
  }
  if (has_interaction) {
    writeLines(interaction_text, conn)
  }
  if (sub_model == "rate") {
    writeLines(timespan_text, conn)
    writeLines(is_dependent_text, conn)
  }
  writeLines(data_text, conn)
  close(conn)

  data_gathered
}
