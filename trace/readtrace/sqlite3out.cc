#include "sqlite3out.h"

#include <stdexcept>

#include "counter.h"
#include "segmentids.h"


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
  auto evtcls = EventclassTable();
  auto evttyp = EventtypeTable(evtcls);
  auto segs = EventTable(evttyp);
  SegmentTable(segs);
}

void SQLiteWriter::LocalityTable() {
  fprintf(stdout, "Adding Locality table..."); fflush(stdout);
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
  fprintf(stdout, "DONE\n"); fflush(stdout);
}

void SQLiteWriter::WorkerTable() {
  fprintf(stdout, "Adding Worker table..."); fflush(stdout);
  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Worker ("
      " id INTEGER, "
      " wid INTEGER, "
      " loc INTEGER, "
      " PRIMARY KEY(id),"
      " FOREIGN KEY (loc) REFERENCES Locality(id));"};
    SQLiteStatement makeit{db_, creation};
    bool more = false;
    do {
      more = makeit.Step();
    } while (more);
  }

  { // INSERT
    std::string addrow{
        "INSERT INTO Worker (id, wid, loc) VALUES (?, ?, ?);"};
    SQLiteStatement addit{db_, addrow};

    Counter worker_id{};

    auto lmap = trace_.worker_map();
    for (auto i = lmap.begin(); i != lmap.end(); ++i) {
      auto wmap = i->second;
      for (auto j = wmap.begin(); j != wmap.end(); ++j) {
        addit.Bind(1, (sqlite3_int64)worker_id.next());
        addit.Bind(2, (sqlite3_int64)j->first);
        addit.Bind(3, (sqlite3_int64)i->first);
        while (addit.Step()) { }
        addit.Reset();
      }
    }
  }
  fprintf(stdout, "DONE\n"); fflush(stdout);
}

void SQLiteWriter::SegmenttypeTable() {
  fprintf(stdout, "Adding Segmenttype table..."); fflush(stdout);

  // Collect data
  std::map<std::string, int> segtyp{};
  trace_.apply([&segtyp](const Event *evt) {
    if (segtyp.find(segment::name(evt->segment_type())) == segtyp.end()) {
      segtyp[segment::name(evt->segment_type())] = evt->segment_type();
    }
  });

  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Segmenttype ("
      " id INTEGER PRIMARY KEY,"
      " name TEXT);"};
    SQLiteStatement makeit{db_, creation};
    while (makeit.Step()) { }
  }

  { // INSERT
    std::string addrow{"INSERT INTO Segmenttype (id, name) VALUES (?, ?);"};
    SQLiteStatement addit{db_, addrow};

    for (auto i = segtyp.begin(); i != segtyp.end(); ++i) {
      addit.Bind(1, (sqlite3_int64)i->second);
      addit.Bind(2, i->first);
      while (addit.Step()) { }
      addit.Reset();
    }
  }

  fprintf(stdout, "DONE\n"); fflush(stdout);
}

std::map<std::string, int> SQLiteWriter::EventclassTable() {
  fprintf(stdout, "Adding Eventclass table..."); fflush(stdout);

  std::map<std::string, int> evtcls{};
  Counter genid{};
  trace_.apply([&evtcls, &genid](const Event *evt) {
    if (evtcls.find(evt->event_class()) == evtcls.end()) {
      evtcls[evt->event_class()] = genid.next();
    }
  });

  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Eventclass ("
      "  id INTEGER,"
      "  class TEXT,"
      "  PRIMARY KEY(id));"};
    SQLiteStatement makeit{db_, creation};
    while (makeit.Step()) { }
  }

  { // INSERT
    std::string addrow{
      "INSERT INTO Eventclass (id, class) VALUES (?, ?);"};
    SQLiteStatement addit{db_, addrow};

    for (auto i = evtcls.begin(); i != evtcls.end(); ++i) {
      addit.Bind(1, (sqlite3_int64)i->second);
      addit.Bind(2, i->first);
      while (addit.Step()) { }
      addit.Reset();
    }
  }

  fprintf(stdout, "DONE\n"); fflush(stdout);

  return evtcls;
}

std::map<std::string, int> SQLiteWriter::EventtypeTable(
                                const std::map<std::string, int> &evtcls) {
  fprintf(stdout, "Adding Eventtype table..."); fflush(stdout);

  std::map<std::string, const Event *> evttyp{};
  trace_.apply([&evttyp](const Event *evt) {
    if (evttyp.find(evt->event_type()) == evttyp.end()) {
      evttyp[evt->event_type()] = evt;
    }
  });

  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Eventtype ("
      " id INTEGER,"
      " cid INTEGER,"
      " name TEXT,"
      " sid INTEGER,"
      " PRIMARY KEY(id),"
      " FOREIGN KEY (sid) REFERENCES Segmenttype(id),"
      " FOREIGN KEY (cid) REFERENCES Eventclass(id));"};
    SQLiteStatement makeit{db_, creation};
    while (makeit.Step()) { }
  }

  std::map<std::string, int> retval{};
  { // INSERT
    std::string addrow{
      "INSERT INTO Eventtype (id, cid, name, sid) VALUES (?, ?, ?, ?);"};
    SQLiteStatement addit{db_, addrow};

    Counter eid{};
    for (auto i = evttyp.begin(); i != evttyp.end(); ++i) {
      retval[i->first] = eid.curr();
      addit.Bind(1, (sqlite3_int64)eid.next());
      addit.Bind(2,
                 (sqlite3_int64)evtcls.find(i->second->event_class())->second);
      addit.Bind(3, i->first);
      addit.Bind(4, (sqlite3_int64)i->second->segment_type());
      while (addit.Step()) { }
      addit.Reset();
    }
  }

  fprintf(stdout, "DONE\n"); fflush(stdout);

  return retval;
}

SegmentCollator SQLiteWriter::EventTable(
      const std::map<std::string, int> &evttyp) {
  fprintf(stdout, "Adding Event table..."); fflush(stdout);

  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Event ("
      " id INTEGER,"
      " tid INTEGER,"
      " timens INTEGER,"
      " PRIMARY KEY(id),"
      " FOREIGN KEY (tid) REFERENCES Eventtype(id));"};
    SQLiteStatement makeit{db_, creation};
    while (makeit.Step()) { }
  }

  SegmentCollator segments{};
  { // INSERT
    std::string addrow{
      "INSERT INTO Event (id, tid, timens) VALUES (?, ?, ?);"};
    SQLiteStatement addit{db_, addrow};

    std::string transtart{"BEGIN TRANSACTION;"};
    SQLiteStatement tstart{db_, transtart};

    std::string tranend{"END TRANSACTION;"};
    SQLiteStatement tend{db_, tranend};

    Counter genid{};
    int tcounter{0};
    trace_.apply([&addit, &genid, &evttyp, &tstart, &tend, &tcounter,
                  &segments] (const Event *evt) {
        if (tcounter == 0) {
          while (tstart.Step()) { }
          tstart.Reset(false);
        }

        if (!segments.collate(evt, genid.curr())) {
          fprintf(stderr, "Event id %lu - Event type '%s'\n", genid.curr(),
                  evt->event_type().c_str());
          throw std::runtime_error("Something with segments.");
        }

        addit.Bind(1, (sqlite3_int64)genid.next());
        addit.Bind(2, (sqlite3_int64)evttyp.find(evt->event_type())->second);
        addit.Bind(3, (sqlite3_int64)evt->stamp());
        while (addit.Step()) { }
        addit.Reset();

        if (++tcounter == 10000) {
          while (tend.Step()) { }
          tend.Reset();
          tcounter = 0;
        }
      }
    );

    if (tcounter > 0 && tcounter < 100000) {
      while (tend.Step()) { }
    }
  }

  fprintf(stdout, "DONE\n"); fflush(stdout);

  return segments;
}

void SQLiteWriter::SegmentTable(const SegmentCollator &segs) {
  fprintf(stdout, "Adding Segment table..."); fflush(stdout);

  { // CREATE TABLE
    std::string creation{
      "CREATE TABLE Segment ("
      " id INTEGER,"
      " sid INTEGER,"
      " start INTEGER,"
      " end INTEGER,"
      " PRIMARY KEY(id),"
      " FOREIGN KEY (start) REFERENCES Event(id),"
      " FOREIGN KEY (end) REFERENCES Event(id));"};
    SQLiteStatement makeit{db_, creation};
    while (makeit.Step()) { }
  }

  { // INSERT
    std::string addrow{
      "INSERT INTO Segment (id, sid, start, end) VALUES (?, ?, ?, ?);"};
    SQLiteStatement addit{db_, addrow};

    std::string transtart{"BEGIN TRANSACTION;"};
    SQLiteStatement tstart{db_, transtart};

    std::string tranend{"END TRANSACTION;"};
    SQLiteStatement tend{db_, tranend};

    int tcounter{0};
    for (size_t i = 0; i < segs.SegmentCount(); ++i) {
      if (tcounter == 0) {
        while (tstart.Step()) { }
        tstart.Reset(false);
      }

      const Segment *curr = segs.GetSegment(i);

      addit.Bind(1, (sqlite3_int64)i);
      addit.Bind(2, (sqlite3_int64)curr->segment());
      addit.Bind(3, (sqlite3_int64)curr->begin());
      addit.Bind(4, (sqlite3_int64)curr->end());
      while (addit.Step()) { }
      addit.Reset();

      if (++tcounter == 10000) {
        while (tend.Step()) { }
        tend.Reset(false);
        tcounter = 0;
      }
    }

    if (tcounter > 0 && tcounter < 100000) {
      while (tend.Step()) { }
    }
  }

  fprintf(stdout, "DONE\n"); fflush(stdout);
}


} // traceutils