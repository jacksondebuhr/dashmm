#include "sqlite3out.h"

#include <stdexcept>

#include "counter.h"


namespace traceutils {


SQLiteStatement::SQLiteStatement(sqlite3 *db, const std::string &sql)
    : db_{db},
      s_{nullptr},
      sql_{sql} {
  // prepare the statement; throw on error
  if (SQLITE_OK != sqlite3_prepare_v2(db_, sql_.c_str(), sql_.length() + 1,
                                      &s_, nullptr)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

SQLiteStatement::~SQLiteStatement() {
  if (s_ != nullptr) {
    sqlite3_finalize(s_);
  }
}

bool SQLiteStatement::Step() {
  int err = sqlite3_step(s_);
  if (err == SQLITE_DONE) {
    return false;
  } else if (err == SQLITE_ROW) {
    return true;
  } else {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
  return false;
}

void SQLiteStatement::Reset(bool clear) {
  if (SQLITE_OK != sqlite3_reset(s_)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
  if (clear && SQLITE_OK != sqlite3_clear_bindings(s_)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

int SQLiteStatement::ParameterCount() {
  return sqlite3_bind_parameter_count(s_);
}

void SQLiteStatement::Bind(int index, double value) {
  if (SQLITE_OK != sqlite3_bind_double(s_, index, value)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

void SQLiteStatement::Bind(int index, sqlite3_int64 value) {
  if (SQLITE_OK != sqlite3_bind_int64(s_, index, value)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

void SQLiteStatement::Bind(int index) {
  if (SQLITE_OK != sqlite3_bind_null(s_, index)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

void SQLiteStatement::Bind(int index, const std::string &value) {
  if (SQLITE_OK != sqlite3_bind_text(s_, index, value.c_str(),
                                       value.length() + 1, SQLITE_TRANSIENT)) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}



SQLiteWriter::SQLiteWriter(const std::string &fname, const Trace &trace)
    : file_name_{fname},
      db_{nullptr},
      trace_{trace} {
  // First we destroy the extant file
  FILE *ifd = fopen(file_name_.c_str(), "r");
  if (ifd != nullptr) {
    throw std::runtime_error("db file already exists. Aborting.");
  }
  //open the connection; throw on error
  int err = sqlite3_open_v2(file_name_.c_str(), &db_,
                            SQLITE_OPEN_READWRITE |SQLITE_OPEN_CREATE,
                            nullptr);
  if (SQLITE_OK != err) {
    throw std::runtime_error(sqlite3_errmsg(db_));
  }
}

SQLiteWriter::~SQLiteWriter() {
  if (db_ != nullptr) {
    sqlite3_close_v2(db_);
  }
}

void SQLiteWriter::write() {
  LocalityTable();
  WorkerTable();
  SegmenttypeTable();
  EventclassTable();
  EventtypeTable();
  EventTable();
  SegmentTable();
}

void SQLiteWriter::LocalityTable() {
  { // CREATE TABLE
    std::string creation{
                  "CREATE TABLE Locality (id INTEGER, PRIMARY KEY(id));"};
    SQLiteStatement makeit{db_, creation};
    bool more = false;
    do {
      more = makeit.Step();
    } while (more);
  }

  { // INSERT
    std::string addrow{"INSERT INTO Locality (id) VALUES (?);"};
    SQLiteStatement addit{db_, addrow};

    // loop over known localities
    for (int i = 0; i <= trace_.max_locality(); ++i) {
      fprintf(stdout, "Got one: %d\n", i);
      if (trace_.has_locality(i)) {
        addit.Bind(1, (sqlite3_int64)i);
        bool more = false;
        do {
          more = addit.Step();
        } while (more);
      }
      addit.Reset();
    }
  }
}

void SQLiteWriter::WorkerTable() {
  //
}

void SQLiteWriter::SegmenttypeTable() {
  //
}

void SQLiteWriter::EventclassTable() {
  //
}

void SQLiteWriter::EventtypeTable() {
  //
}

void SQLiteWriter::EventTable() {
  //
}

void SQLiteWriter::SegmentTable() {
  //
}


} // traceutils